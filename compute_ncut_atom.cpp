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

#include "compute_ncut_atom.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNCutAtom::ComputeNCutAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute ncut/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  double cutoff = force->numeric(FLERR,arg[3]);
  if (cutoff < 0.0) error->all(FLERR,"Illegal compute ncut/atom command");
  cutsq = cutoff*cutoff;

  nmax = 0;
  ncut = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeNCutAtom::~ComputeNCutAtom()
{
  memory->destroy(ncut);
}

/* ---------------------------------------------------------------------- */

void ComputeNCutAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute ncut/atom requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute ncut/atom cutoff is longer than pairwise cutoff");

  // cannot use neighbor->cutneighmax b/c neighbor has not yet been init

  if (2.0*sqrt(cutsq) > force->pair->cutforce + neighbor->skin &&
      comm->me == 0)
    error->warning(FLERR,"Compute ncut/atom cutoff may be too large to find "
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

void ComputeNCutAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeNCutAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;

  invoked_peratom = update->ntimestep;

  // grow arrays if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(ncut);
    nmax = atom->nmax;

    memory->create(ncut,nmax,"ncut/atom:ncut");
    vector_atom = ncut;
  }

  // invoke full neighbor list (will copy or build if necessary)

  //XXX Again, use this flag so that we can change_box, run 0, and have updates
  //XXX                            V
  neighbor->build_one(list->index, 1);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
	if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

	ncut[i] = 0;
    for (jj = 0; jj < jnum; ++jj) {
      j = jlist[jj];
	  if (!(mask[j] & groupbit)) continue;
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) ++ncut[i];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeNCutAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
