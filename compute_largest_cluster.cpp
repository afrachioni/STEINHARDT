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
#include "compute_largest_cluster.h"
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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeLargestCluster::ComputeLargestCluster(LAMMPS *lmp, int narg, char **arg) :
	Compute(lmp, narg, arg)
{
	if (narg != 4) error->all(FLERR, "Illegal compute largest cluster command");

	scalar_flag = 1;

	int icompute = modify->find_compute(arg[3]);
	source_compute = modify->compute[icompute];
	if (icompute < 0)
		error->all(FLERR, "Compute ID for compute largest cluster does not exist");
	if (!source_compute->peratom_flag)
		error->all(FLERR, "Compute specified to compute largest_cluster is not per atom");
	if (source_compute->size_peratom_cols)
		error->all(FLERR, "Compute specified to compute largest_cluster does not compute a per-atom vector");
}

/* ---------------------------------------------------------------------- */

ComputeLargestCluster::~ComputeLargestCluster()
{
}

/* ---------------------------------------------------------------------- */

void ComputeLargestCluster::init()
{
}

/* ---------------------------------------------------------------------- */

//void ComputeLargestCluster::init_list(int id, NeighList *ptr)
//{
	//list = ptr;
//}

/* ---------------------------------------------------------------------- */

double ComputeLargestCluster::compute_scalar()
{
	invoked_scalar = update->ntimestep;

	// This is a very naive, not scalable implementation so that I can get up and running
	// watch out!
	if (!source_compute->invoked_flag) // TODO do we need timestep logic?
		source_compute->compute_peratom();

	double *source = source_compute->vector_atom;
	int nlocal = atom -> nlocal;
	int *mask = atom -> mask;

	double global_clusters[atom->natoms];
	double cluster_hist[atom->natoms];
	for (int i = 0; i < atom->natoms; ++i)
		cluster_hist[i] = 0;

	int nlocal_list[comm->nprocs];
	int offsets[comm->nprocs];
	MPI_Gather (&nlocal, 1, MPI_INT, nlocal_list, 1, MPI_INT, 0, world);

	offsets[0] = 0;
	for (int i = 1; i < comm->nprocs; ++i)
		offsets[i] = offsets[i - 1] + nlocal_list[i - 1];

	MPI_Gatherv (source, nlocal, MPI_DOUBLE, global_clusters, nlocal_list, offsets, MPI_DOUBLE, 0, world);

	double max = -1;
	if (comm->me == 0) {
		for (int i = 0; i < atom->natoms; ++i) {
			if (global_clusters[i] != 0)
				++(cluster_hist[(int)global_clusters[i] - 1]); // one based indicies
		}

		for (int i = 0; i < atom->natoms; ++i)
			if (cluster_hist[i] > max)
				max = cluster_hist[i];
	}
	MPI_Bcast (&max, 1, MPI_DOUBLE, 0, world);

	scalar = max;
	return scalar;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
   ------------------------------------------------------------------------- */

double ComputeLargestCluster::memory_usage()
{
	return 0;
}
