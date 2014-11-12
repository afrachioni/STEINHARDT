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

#ifdef COMPUTE_CLASS

ComputeStyle(qdotq/atom,ComputeQDotQ)

#else

#ifndef LMP_COMPUTE_QDOTQ_H
#define LMP_COMPUTE_QDOTQ_H

#include "compute.h"

namespace LAMMPS_NS {
	class ComputeQDotQ : public Compute {
		public:
			ComputeQDotQ(class LAMMPS *, int, char **);
			~ComputeQDotQ();
			void init();
			void init_list(int, class NeighList *);
			void compute_peratom();
			double memory_usage();
			private:
			Compute *source_compute;

			int nmax;
			int nghost_max;
			double cutsq;
			class NeighList *list;
			double **qlm_buffer;
			int **nearest;
			int *nnearest;
			double *pattern;

			int l;
			double threshold;
	};
}

#endif
#endif
