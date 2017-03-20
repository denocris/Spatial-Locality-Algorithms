#ifndef __TOOLS__
#define __TOOLS__

#include <string>
#include <iostream>
#include <cmath>

typedef struct Hilbert{
  int npoints;
  double * xpos;
  double * ypos;
  double **  quad_pos;
  double **  corner;
  int reclevel, recdepth, side;

} hcurve;

void hilbert_lattice_init(Hilbert & hilb,
                          int,
                          int);

void hilbert_dealloc(Hilbert & hilb);

void make_hilbert_grid(Hilbert & hilb,
                       double * xgrid,
                       double * ygrid,
                       int n);

#endif
