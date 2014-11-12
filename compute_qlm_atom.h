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

ComputeStyle(qlm/atom,ComputeQLMAtom)

#else

#ifndef LMP_COMPUTE_QLM_ATOM_H
#define LMP_COMPUTE_QLM_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {
	class ComputeQLMAtom : public Compute {
		public:
			ComputeQLMAtom(class LAMMPS *, int, char **);
			~ComputeQLMAtom();
			void init();
			void init_list(int, class NeighList *);
			void compute_peratom();
			double memory_usage();

		private:
			int l;
			int nghost_max;
			double cutsq;
			class NeighList *list;
			int **nearest;
			int *nnearest;
			double **qlm_buffer;
	};
}

#endif
#endif
