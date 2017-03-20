#ifndef __TOOLS__
#define __TOOLS__

#include <string>
#include <iostream>
#include <cmath>
#include <bitset>

typedef std::bitset<64> bittype; // Alias

typedef struct Morton{
  double * xpos;
  double * ypos;
  double ** quad;
  double ** corner;
  int reclevel, recdepth, size, npoints;
  int * key_list;
} mcurve;

bool mybitf(bittype & n, int pos);

void mybitset(bittype & n, int pos);

void morton_lattice_init(Morton & mort, int scan, int size);

void morton_dealloc(Morton & mort);

void make_morton_grid(Morton & mort,
                      double * xgrid,
                      double * ygrid,
                      int n);

void make_morton_key(Morton* mort);

#endif
