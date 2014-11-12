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
   Contributing author: Anthony Frachioni
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "compute_boop_atom.h"
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

#define MAXNEAR 20

enum{NCOMMON,NBOND,MAXBOND,MINBOND};

/* ---------------------------------------------------------------------- */

ComputeBOOPAtom::ComputeBOOPAtom(LAMMPS *lmp, int narg, char **arg) :
	Compute(lmp, narg, arg)
{
	if (narg != 5) error->all(FLERR, "Illegal compute BOOP command");

	peratom_flag = 1;
	size_peratom_cols = 0;

	l = atoi(arg[3]);//TODO check behavior if arg[3] not an integer
	//TODO check l>8
	if (l < 1 || l > 8) error->all(FLERR, "Illegal compute BOOP command");

	double cutoff = atof(arg[4]);
	if (cutoff < 0.0) error->all(FLERR, "Illegal compute BOOP command");
	cutsq = cutoff*cutoff;

	nmax = 0;
	nearest = NULL;
	nnearest = NULL;
	pattern = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBOOPAtom::~ComputeBOOPAtom()
{
	memory->destroy(nearest);
	memory->destroy(nnearest);
	memory->destroy(pattern);
}

/* ---------------------------------------------------------------------- */

void ComputeBOOPAtom::init()
{
	if (force->pair == NULL) 
		error->all(FLERR, "Compute boop requires a pair style be defined");
	if (sqrt(cutsq) > force->pair->cutforce) 
		error->all(FLERR, "Compute boop cutoff is longer than pairwise cutoff");

	// cannot use neighbor->cutneighmax b/c neighbor has not yet been init

	if (2.0*sqrt(cutsq) > force->pair->cutforce + neighbor->skin && comm->me == 0)
		error->warning(FLERR, "Compute boop cutoff may be too large to find " "ghost atom neighbors");

	// need an occasional full neighbor list
	int irequest = neighbor->request((void *) this);
	neighbor->requests[irequest]->pair = 0;
	neighbor->requests[irequest]->compute = 1;
	neighbor->requests[irequest]->half = 0;
	neighbor->requests[irequest]->full = 1;
	neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeBOOPAtom::init_list(int id, NeighList *ptr)
{
	list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeBOOPAtom::compute_peratom()
{
	int i,j,ii,jj,n,inum,jnum;
	int *ilist,*jlist,*numneigh,**firstneigh;
	double xtmp,ytmp,ztmp,delx,dely,delz,rsq;

	invoked_peratom = update->ntimestep;

	// grow arrays if necessary
	if (atom->nlocal > nmax) {
		memory->destroy(nearest);
		memory->destroy(nnearest);
		memory->destroy(pattern);
		nmax = atom->nmax;

		memory->create(nearest,nmax,MAXNEAR,"boop:nearest");
		memory->create(nnearest,nmax,"boop:nnearest");
		memory->create(pattern,nmax,"boop:boop_pattern");
		vector_atom = pattern;
		//TODO: rename pattern something useful
	}

	// invoke full neighbor list (will copy or build if necessary)
	neighbor->build_one(list);

	inum = list->inum;
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;

	// find the neigbours of each atom within cutoff using full neighbor list
	// nearest[] = atom indices of nearest neighbors, up to MAXNEAR
	// do this for all atoms, not just compute group
	// since CNA calculation requires neighbors of neighbors

	double **x = atom->x;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;

	int nerror = 0;
	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		if (!(mask[i] & groupbit)) continue;
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

	int nerrorall;
	MPI_Allreduce(&nerror,&nerrorall,1,MPI_INT,MPI_SUM,world);
	if (nerrorall && comm->me == 0) {
		char str[128];
		sprintf(str,"Too many neighbors in BOOP for %d atoms",nerrorall);
		error->warning(FLERR, str,0);
	}

	//Compute BOOP per atom
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int current_index;
	int *his_neighbors;
	double d_x, d_y, d_z, theta, phi;
	double current_x, current_y, current_z;
	double ql;
	std::complex<double> qlm;
	double pi = boost::math::constants::pi<double>();
	for (int iii = 0; iii < inum; iii++) { 
		current_index = ilist[iii];
		if (!(mask[current_index] & groupbit)) continue;
		his_neighbors = nearest[current_index];
		current_x = x[current_index][0];
		current_y = x[current_index][1];
		current_z = x[current_index][2];
		ql = 0;
		for (int m = -l; m < l + 1; m++) {
			qlm = std::complex<double> (0, 0);
			for (int j = 0; j < nnearest[current_index]; j++) {
				d_x = current_x - x[his_neighbors[j]][0];
				d_y = current_y - x[his_neighbors[j]][1];
				d_z = current_z - x[his_neighbors[j]][2];
				theta = acos (d_z / sqrt(d_x*d_x + d_y*d_y + d_z*d_z));
				phi = atan2 (d_y, d_x);
				qlm += boost::math::spherical_harmonic(l, m, theta, phi);
			}
			qlm /= nnearest[current_index];
			ql += norm (qlm);
		}
		ql *= 4 * pi / (2 * l + 1);
		pattern[current_index] = sqrt(ql);
	}
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   ------------------------------------------------------------------------- */

double ComputeBOOPAtom::memory_usage()
{
	double bytes = nmax * sizeof(int);
	bytes += nmax * MAXNEAR * sizeof(int);
	bytes += nmax * sizeof(double);
	return bytes;
}
