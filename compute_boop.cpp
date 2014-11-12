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
#include "compute_boop.h"
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

//#define MAXNEAR 16
#define MAXNEAR 100 //TODO I don't see the point in keeping this around
#define MAXCOMMON 8

enum{NCOMMON,NBOND,MAXBOND,MINBOND};

/* ---------------------------------------------------------------------- */

ComputeBOOP::ComputeBOOP(LAMMPS *lmp, int narg, char **arg) :
	Compute(lmp, narg, arg)
{
	if (narg != 5) error->all(FLERR, "Illegal compute BOOP command");

	scalar_flag = 1;

	l = atoi(arg[3]);//TODO check behavior if arg[3] not an integer
	if (l < 1 || l > 8) error->all(FLERR, "Illegal compute BOOP command");
	//TODO see if and when overflows start happening in l

	double cutoff = atof(arg[4]);
	if (cutoff < 0.0) error->all(FLERR, "Illegal compute BOOP command");
	cutsq = cutoff*cutoff;

	nmax = 0;
	nearest = NULL;
	nnearest = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBOOP::~ComputeBOOP()
{
	memory->destroy(nearest);
	memory->destroy(nnearest);
}

/* ---------------------------------------------------------------------- */

void ComputeBOOP::init()
{
	if (force->pair == NULL) 
		error->all(FLERR, "Compute boop requires a pair style be defined");
	if (sqrt(cutsq) > force->pair->cutforce) 
		error->all(FLERR, "Compute boop cutoff is longer than pairwise cutoff");

	// cannot use neighbor->cutneighmax b/c neighbor has not yet been init
	if (2.0*sqrt(cutsq) > force->pair->cutforce + neighbor->skin && comm->me == 0)
		error->warning(FLERR, "Compute boop cutoff may be too large to find "
				"ghost atom neighbors");

	// need an occasional full neighbor list
	int irequest = neighbor->request((void *) this);
	neighbor->requests[irequest]->pair = 0;
	neighbor->requests[irequest]->compute = 1;
	neighbor->requests[irequest]->half = 0;
	neighbor->requests[irequest]->full = 1;
	neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeBOOP::init_list(int id, NeighList *ptr)
{
	list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeBOOP::compute_scalar()
{
	int i,j,ii,jj,n,inum,jnum;
	int *ilist,*jlist,*numneigh,**firstneigh;
	double xtmp,ytmp,ztmp,delx,dely,delz,rsq;

	invoked_scalar = update->ntimestep;

	// grow arrays if necessary
	if (atom->nlocal > nmax) {
		memory->destroy(nearest);
		memory->destroy(nnearest);
		nmax = atom->nmax;

		memory->create(nearest,nmax,MAXNEAR,"boop:nearest");
		memory->create(nnearest,nmax,"boop:nnearest");
	}

	// invoke full neighbor list (will copy or build if necessary)
	// XXX This flag can be set to zero for ordinary use, and maybe should
	// XXX exposed to the user, but  |  is necessary for hybrid MC NPT ensemble
	// XXX                           V
	neighbor->build_one(list, 1);

	inum = list->inum;
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;


	// find the neigbours of each atom within cutoff using full neighbor list
	double **x = atom->x;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;

	int nerror = 0;

	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		// is this proper?
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

	//Compute BOOP
	int current_index, num_neighbors;
	int zero_neighbors_count = 0;
	int *his_neighbors;
	double d_x, d_y, d_z, theta, phi;
	double current_x, current_y, current_z;
	std::complex<double> qlm, qlm_local_sum;
	double qlm_global_sum_real, qlm_global_sum_imag;
	double qlm_local_sum_real, qlm_local_sum_imag;
	double pi = boost::math::constants::pi<double>();
	double factor = 4 * pi / (2 * l + 1);
	double m_sum = 0;
	for (int m = -l; m < l + 1; m++) {
		qlm_local_sum = std::complex<double> (0, 0);
		for (int iii = 0; iii < inum; iii++) {
			current_index = ilist[iii];
			// XXX note: terms are added for neighbors not in the compute group
			// but I think this is OK
			if (!(mask[current_index] & groupbit)) continue;
			num_neighbors = nnearest[current_index];
			if (num_neighbors == 0) {
				zero_neighbors_count++;
				continue;
			}
			his_neighbors = nearest[current_index];
			current_x = x[current_index][0];
			current_y = x[current_index][1];
			current_z = x[current_index][2];
			qlm = std::complex<double> (0, 0);
			for (int neigh_index = 0; neigh_index < num_neighbors; neigh_index++) {
				d_x = current_x - x[his_neighbors[neigh_index]][0];
				d_y = current_y - x[his_neighbors[neigh_index]][1];
				d_z = current_z - x[his_neighbors[neigh_index]][2];
				theta = acos (d_z / sqrt( d_x * d_x + d_y * d_y + d_z * d_z));	
				phi = atan2 (d_y, d_x);
				qlm += boost::math::spherical_harmonic(l, m, theta, phi);
			}
			qlm /= num_neighbors;
			qlm_local_sum += qlm;
		}
		qlm_local_sum_real = real(qlm_local_sum);
		qlm_local_sum_imag = imag(qlm_local_sum);
		MPI_Allreduce ( &qlm_local_sum_real, &qlm_global_sum_real, 1, MPI_DOUBLE, MPI_SUM, world);
		MPI_Allreduce ( &qlm_local_sum_imag, &qlm_global_sum_imag, 1, MPI_DOUBLE, MPI_SUM, world);
		m_sum += (qlm_global_sum_real*qlm_global_sum_real + qlm_global_sum_imag*qlm_global_sum_imag);
	}
	scalar = sqrt (factor * m_sum) / (atom->natoms - zero_neighbors_count);
	return scalar;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   ------------------------------------------------------------------------- */

double ComputeBOOP::memory_usage()
{
	double bytes = nmax * sizeof(int);
	bytes += nmax * MAXNEAR * sizeof(int);
	bytes += nmax * sizeof(double);
	return bytes;
}
