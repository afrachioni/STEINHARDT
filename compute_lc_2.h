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

ComputeStyle(lc2,ComputeLC2)

#else

#ifndef LMP_COMPUTE_LC2_H
#define LMP_COMPUTE_LC2_H

#include "compute.h"

namespace LAMMPS_NS {
	class ComputeLC2 : public Compute {
		public:
			ComputeLC2(class LAMMPS *, int, char **);
			~ComputeLC2();
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
