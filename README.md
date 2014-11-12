This is a package for LAMMPS, which includes a few compute classes for calculation of Steinhardt parameters and various related quantities.  This package was last successfully built against the **5 Nov 2014** version of LAMMPS.  Compute-specific information follows:

---


####`compute boop` and `compute boop/atom`:
Implementation of the Steinhardt Q parameters, for an entire group and per-atom.

Stable, but depends on external BOOST for spherical harmonics.

Usage:
```
compute ID group-ID boop l cutoff
compute ID group-ID boop/atom l cutoff
```

---


####`compute ncut/atom`:
Number of atoms within some cutoff radius, per-atom.  This surprisingly seems missing from LAMMPS.

Usage:
```
compute ID group-ID ncut/atom cutoff
```


---


####`compute largest_cluster`:
Not stable or fast.  Calculates the number of atoms in the largest cluster given the ID of a cluster/atom compute.

Usage:
```
compute ID largest_cluster cluster_atom-ID
```
where `cluster_atom-ID` is the ID of a compute of style cluster/atom.


---


####`compute qlm/atom`:
Might work.  Normalized *q<sup>(i)</sup><sub>lm</sub>*, see Wolde 96 [15].  Output is a per-atom array with 4 * l + 2 columns (2 * l + 1 complex numbers).

Usage:
```
compute ID group-ID qlm/atom l cutoff
```


---


####`compute qdotq/atom`:
*l* hardcoded to 6.  Number of neighbors with ||*q<sup>(i)</sup><sub>l</sub>*  ·  *q<sup>(j)</sup><sub>l</sub>*|| above some threshold.  See Wolde 96 [16].

Usage:
```
compute ID group-ID qdotq/atom qlm-ID cutoff threshold
```
where `qlm-ID` is the ID of a compute of style qlm/atom.


---


####`compute cluster_qdotq/atom`:
*l* hardcoded to 6.  As with `compute cluster/atom`, assigns each atom a cluster ID.  Atoms belong to the same cluster if they are less than `cutoff` apart, and the dot product *q<sup>(i)</sup><sub>l</sub>*  ·  *q<sup>(j)</sup><sub>l</sub>* exceeds `threshold`.  See Wolde 96.

Usage:
```
compute ID group-ID cluster_qdotq/atom cutoff qdotq-ID threshold
```
where `qlm-ID` is the ID of a compute of style qlm/atom.
