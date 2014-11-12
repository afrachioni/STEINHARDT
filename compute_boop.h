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

#ifdef COMPUTE_CLASS

ComputeStyle(boop,ComputeBOOP)

#else

#ifndef LMP_COMPUTE_BOOP_H
#define LMP_COMPUTE_BOOP_H

#include "compute.h"

namespace LAMMPS_NS {
	class ComputeBOOP : public Compute {
		public:
			ComputeBOOP(class LAMMPS *, int, char **);
			~ComputeBOOP();
			void init();
			void init_list(int, class NeighList *);
			double compute_scalar();
			double memory_usage();

		private:
			int l;
			int nmax;
			double cutsq;
			class NeighList *list;
			int **nearest;
			int *nnearest;
	};
}

#endif
#endif
