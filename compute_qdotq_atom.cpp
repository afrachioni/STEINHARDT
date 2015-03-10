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
#include "compute_qdotq_atom.h"
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

/* ---------------------------------------------------------------------- */

ComputeQDotQ::ComputeQDotQ(LAMMPS *lmp, int narg, char **arg) :
	Compute(lmp, narg, arg)
{
	// Usage: compute ID group-ID qdotq/atom qlm-ID cutoff threshold

	if (narg != 6) error->all(FLERR, "Illegal compute largest cluster command");

	peratom_flag = 1;
	size_peratom_cols = 0;

	l = 6; // XXX

	/*
	int icompute = modify->find_compute(arg[3]);
	source_compute = modify->compute[icompute];
	if (icompute < 0)
		error->all(FLERR, "Compute ID for compute qdotq does not exist");
	if (strcmp(source_compute->style, "qlm/atom"))
		error->all(FLERR, "Compute specified to compute qdotq is not of style qlm/atom");
*/

	double cutoff = atof(arg[4]);
	if (cutoff < 0) error->all(FLERR, "Illegal compute qdotq/atom");
	cutsq = cutoff*cutoff;

	threshold = atof(arg[5]);
	if (threshold < 0) error->all(FLERR, "Illegal compute qdotq/atom");

	nmax = 0;
	nghost_max = 0;
	qlm_buffer = NULL;
	nearest = NULL;
	nnearest = NULL;
	pattern = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeQDotQ::~ComputeQDotQ()
{
	memory->destroy(qlm_buffer);
	memory->destroy(nearest);
	memory->destroy(nnearest);
	memory->destroy(pattern);
}

/* ---------------------------------------------------------------------- */

void ComputeQDotQ::init()
{
	if (force->pair == NULL)
		error->all(FLERR, "Compute qdotq requires a pair style be defined");

	if (sqrt(cutsq) > force->pair->cutforce)
		error->all(FLERR, "Compute qdotq/atom cutoff is longer than pairwise cutoff");

	if (2*sqrt(cutsq) > force->pair->cutforce + neighbor->skin && comm->me == 0)
		error->warning(FLERR, "Compute qdotq/atom cutoff may be too large to "
				"find ghost atom neighbors");

	// Occasional full neighbor list, with ghosts
	int irequest = neighbor->request((void *) this);
	neighbor->requests[irequest]->pair = 0;
	neighbor->requests[irequest]->compute = 1;
	neighbor->requests[irequest]->half = 0;
	neighbor->requests[irequest]->full = 1;
	neighbor->requests[irequest]->occasional = 1; // This was zero for some reason
	neighbor->requests[irequest]->ghost = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeQDotQ::init_list(int id, NeighList *ptr)
{
	list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeQDotQ::compute_peratom()
{
	int i, j, ii, jj, n, inum, jnum, gnum;
	int *ilist, *jlist, *numneigh, **firstneigh;
	double xtmp, ytmp, ztmp, delx, dely, delz, rsq;

	invoked_peratom = update->ntimestep;

	//if (!source_compute->invoked_flag) // TODO do we need timestep logic?
		//source_compute->compute_peratom();

	// grow arrays if necessary
	if (atom->nlocal + atom->nghost > nghost_max) {
		memory->destroy(nearest);
		memory->destroy(nnearest);
		memory->destroy(qlm_buffer);
		nghost_max = atom->nlocal + atom->nghost;

		memory->create(nearest, nghost_max, MAXNEAR, "qdotq:nearest");
		memory->create(nnearest, nghost_max, "qdotq:nnearest");
		memory->create(qlm_buffer, nghost_max, 4 * l + 2, "qdotq:qlm_buffer");
	}
	if (atom->nlocal > nmax) {
		memory->destroy(pattern);
		nmax = atom->nmax;
		memory->create(pattern, nmax, "qdotq:pattern");
		vector_atom = pattern;
	}

	neighbor->build_one(list);

	inum = list->inum;
	gnum = list->gnum;
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;

	//find neighbors

	double **x = atom -> x;
	int *mask = atom->mask;
	int nlocal = atom ->nlocal;

	int nerror = 0;

	for (ii = 0; ii < inum + gnum; ++ii) {
		i = ilist[ii];
		if (!(mask[i] & groupbit)) continue;
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];
		jlist = firstneigh[i];
		jnum = numneigh[i];

		n = 0;
		for (jj = 0; jj < jnum; ++jj) {
			j = jlist[jj];
			j &= NEIGHMASK;

			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			rsq = delx * delx + dely * dely + delz * delz;
			if (rsq < cutsq) {
				if (n < MAXNEAR) nearest[i][n++] = j;
				else {
					++nerror;
					break;
				}
			}
		}
			nnearest[i] = n;
	}
	int nerrorall;
	MPI_Allreduce (&nerror, &nerrorall, 1, MPI_INT, MPI_SUM, world);
	if (nerrorall && comm->me == 0) {
		char str[128];
		sprintf (str, "Too many neighbors in qdotq for %d atoms", nerrorall);
		error->warning (FLERR, str, 0);
	}

	// Compute qlm per atom, including ghosts (boo!)
	int current_index;
	double d_x, d_y, d_z, theta, phi;
	double current_x, current_y, current_z;
	double norm;
	std::complex<double> qlm;
	double pi = boost::math::constants::pi<double>();
	for (int iii = 0; iii < inum + gnum; ++iii) {
		current_index = ilist[iii];
		current_x = x[current_index][0];
		current_y = x[current_index][1];
		current_z = x[current_index][2];

		for (int m = -l; m < l + 1; ++m) {
			qlm = std::complex<double> (0, 0);
			for (int jjj = 0; jjj < nnearest[current_index]; ++jjj) {
				j = nearest[current_index][jjj];
				j &= NEIGHMASK;
				d_x = current_x - x[j][0];
				d_y = current_y - x[j][1];
				d_z = current_z - x[j][2];

				theta = acos (d_z / sqrt(d_x*d_x + d_y*d_y + d_z*d_z));
				phi = atan2 (d_y, d_x);
				qlm += boost::math::spherical_harmonic(l, m, theta, phi);
			}
			if (1 || nnearest[current_index] != 0) //XXX
				qlm /= nnearest[current_index];
			qlm_buffer[current_index][2 * (l + m)] = real(qlm);
			qlm_buffer[current_index][2 * (l + m) + 1] = imag(qlm);
		}
		norm = 0;
		for (int m = -l; m < l + 1; ++m) // normalize! (maybe this would be clearer with different loop)
			norm += qlm_buffer[current_index][2*(l+m)] * qlm_buffer[current_index][2*(l+m)] + \
				qlm_buffer[current_index][2*(l+m)+1] * qlm_buffer[current_index][2*(l+m)+1];
		norm = sqrt(norm);
		for (int k = 0; k < 4 * l + 2; ++k)
			if (1 || qlm_buffer[current_index][k] != 0)
				qlm_buffer[current_index][k] /= norm;
	}


	// Execute dot product
	for (int iii = 0; iii < inum + gnum; ++iii) {
		int count = 0;
		current_index = ilist[iii];
		if (!(mask[current_index] & groupbit)) continue;
		
		double *my_qlm = qlm_buffer[current_index];

		for (int j = 0; j < nnearest[current_index]; ++j) {
			double dot = 0;
			int kk = nearest[current_index][j];
			double *his_qlm = qlm_buffer[nearest[current_index][j]];
			double my_real, his_real, my_imaj, his_imaj;
			for (int m = -l; m < l + 1; ++m) {
				my_real = my_qlm[2 * (m + l)];
				his_real = his_qlm[2 * (m +l)];
				my_imaj = my_qlm[2 * (m + l) + 1];
				his_imaj = his_qlm[2 * (m + l) + 1];

				dot += my_real*his_real + my_imaj*his_imaj;
			}
			if (dot > threshold)
				++count;
		}
		pattern[current_index] = count;
	}
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   ------------------------------------------------------------------------- */

double ComputeQDotQ::memory_usage()
{
	return 0;
}
