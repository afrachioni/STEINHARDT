/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Anthony Frachioni (Binghamton University)
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "compute_qlm_atom.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math.h"

#include<boost/math/special_functions/spherical_harmonic.hpp>

using namespace LAMMPS_NS;

#define MAXNEAR 24
#define MAXCOMMON 8

enum{NCOMMON,NBOND,MAXBOND,MINBOND};

/* ---------------------------------------------------------------------- */

ComputeQLMAtom::ComputeQLMAtom(LAMMPS *lmp, int narg, char **arg) :
	Compute(lmp, narg, arg)
{
	if (narg != 5) error->all(FLERR, "Illegal compute BOOP command");

	l = atoi(arg[3]);//TODO check behavior if arg[3] not an integer
	//TODO check l>8
	if (l < 1 || l > 8) error->all(FLERR, "Illegal compute BOOP command");

	peratom_flag = 1;
	size_peratom_cols = 4 * l + 2;

	double cutoff = atof(arg[4]);
	if (cutoff < 0.0) error->all(FLERR, "Illegal compute BOOP command");
	cutsq = cutoff*cutoff;

	nghost_max = 0;
	nghost_max = 0;
	nearest = NULL;
	nnearest = NULL;
	qlm_buffer = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeQLMAtom::~ComputeQLMAtom()
{
	memory->destroy(nearest);
	memory->destroy(nnearest);
	memory->destroy(qlm_buffer);
}

/* ---------------------------------------------------------------------- */

void ComputeQLMAtom::init()
{
	if (force->pair == NULL) 
		error->all(FLERR, "Compute qlm_atom requires a pair style be defined");
	if (sqrt(cutsq) > force->pair->cutforce) 
		error->all(FLERR, "Compute qlm_atom cutoff is longer than pairwise cutoff");

	// cannot use neighbor->cutneighmax b/c neighbor has not yet been init

	if (2.0*sqrt(cutsq) > force->pair->cutforce + neighbor->skin && comm->me == 0)
		error->warning(FLERR, "Compute qlm_atom cutoff may be too large to find " "ghost atom neighbors");

	int count = 0;
	for (int i = 0; i < modify->ncompute; i++)
		if (strcmp(modify->compute[i]->style,"qlm_atom") == 0) count++;
	if (count > 1 && comm->me == 0)
		error->warning(FLERR, "More than one compute qlm_atom defined");
	//TODO multiple qlm_atoms is probably ok

	// need an occasional full neighbor list

	int irequest = neighbor->request((void *) this);
	neighbor->requests[irequest]->pair = 0;
	neighbor->requests[irequest]->compute = 1;
	neighbor->requests[irequest]->half = 0;
	neighbor->requests[irequest]->full = 1;
	neighbor->requests[irequest]->occasional = 1;
	// XXX with ghosts
	neighbor->requests[irequest]->ghost = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeQLMAtom::init_list(int id, NeighList *ptr)
{
	list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeQLMAtom::compute_peratom()
{
	int i,j,ii,jj,n,inum,jnum,gnum;
	int *ilist,*jlist,*numneigh,**firstneigh;
	double xtmp,ytmp,ztmp,delx,dely,delz,rsq;

	invoked_peratom = update->ntimestep;

	// grow arrays if necessary
	if (atom->nlocal + atom->nghost > nghost_max) {
		memory->destroy(nearest);
		memory->destroy(nnearest);
		memory->destroy(qlm_buffer);
		nghost_max = atom->nlocal + atom->nghost;

		memory->create(nearest,nghost_max,MAXNEAR,"qlm_atom:nearest");
		memory->create(nnearest,nghost_max,"qlm_atom:nnearest");
		memory->create(qlm_buffer, nghost_max, 4 * l + 2, "qlm_atom:qlm_buffer");
		array_atom = qlm_buffer;
	}

	// invoke full neighbor list (will copy or build if necessary)
	neighbor->build_one(list->index);

	inum = list->inum;
	gnum = list->gnum;
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;

	// find the neigbours of each atom within cutoff using full neighbor list
	// nearest[] = atom indices of nearest neighbors, up to MAXNEAR
	// do this for all atoms, not just compute group
	// since CNA calculation requires neighbors of neighbors
	//TODO do this only for compute group

	double **x = atom->x;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;

	int nerror = 0;
	for (ii = 0; ii < inum + gnum; ++ii) {
		i = ilist[ii];
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];
		jlist = firstneigh[i];
		jnum = numneigh[i];

		n = 0;
		for (jj = 0; jj < jnum; jj++) {
			j = jlist[jj];
			j &= NEIGHMASK;

			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			rsq = delx*delx + dely*dely + delz*delz;
			if (rsq < cutsq) {
				if (n < MAXNEAR) nearest[i][n++] = j;
				else {
					nerror++;
					break;
				}
			}
		}
		nnearest[i] = n;
	}

	// warning message
	//TODO handle no neighbors case, maybe with warning
	int nerrorall;
	MPI_Allreduce(&nerror,&nerrorall,1,MPI_INT,MPI_SUM,world);
	if (nerrorall && comm->me == 0) {
		char str[128];
		sprintf(str,"Too many neighbors in qlm for %d atoms",nerrorall);
		error->warning(FLERR, str,0);
	}

	//Compute BOOP per atom, including ghosts
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int current_index;
	double d_x, d_y, d_z, theta, phi;
	double current_x, current_y, current_z;
	double norm;
	std::complex<double> qlm;
	double pi = boost::math::constants::pi<double>();
	for (int iii = 0; iii < inum + gnum; iii++) { 
		current_index = ilist[iii];
		current_x = x[current_index][0];
		current_y = x[current_index][1];
		current_z = x[current_index][2];
		// zero buffer
		for (int k = 0; k < 4 * l + 2; ++k)
			qlm_buffer[current_index][k] = 0;
		for (int m = -l; m < l + 1; m++) {
			qlm = std::complex<double> (0, 0);
			for (int jj = 0; jj < nnearest[current_index]; ++jj) {
				j = nearest[current_index][jj];
				j &= NEIGHMASK;
				d_x = current_x - x[j][0];
				d_y = current_y - x[j][1];
				d_z = current_z - x[j][2];

				theta = acos (d_z / sqrt(d_x*d_x + d_y*d_y + d_z*d_z));
				phi = atan2 (d_y, d_x);
				qlm += boost::math::spherical_harmonic(l, m, theta, phi);
			}
			qlm /= nnearest[current_index];
			qlm_buffer[current_index][2 * (l + m)] = real(qlm);
			qlm_buffer[current_index][2 * (l + m) + 1] = imag(qlm);
		}
		norm = 0;
		for (int m = -l; m < l + 1; ++m)
			norm += qlm_buffer[current_index][2*(l+m)] * qlm_buffer[current_index][2*(l+m)] + \
				qlm_buffer[current_index][2*(l+m)+1] * qlm_buffer[current_index][2*(l+m)+1];
		norm = sqrt(norm);
		for (int k = 0; k < 4 * l + 2; ++k)
			if (qlm_buffer[current_index][k] != 0)
					qlm_buffer[current_index][k] /= norm;
	}
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   ------------------------------------------------------------------------- */

double ComputeQLMAtom::memory_usage()
{
	double bytes = nghost_max * sizeof(int);
	bytes += nghost_max * MAXNEAR * sizeof(int);
	bytes += nghost_max * (4 * l + 2) * sizeof(double);
	return bytes;
}
