#include "tools.h"


int main(){

  int i, j, niterations;
  double * xgrid, * ygrid;
  FILE * out;
  struct Morton mcurve;

  xgrid = new double [4];
  ygrid = new double [4];

  niterations = 4;

  morton_lattice_init(mcurve, niterations, 16);

/* ------------------- Initialization ---------------- */
  for(j = 0; j < 4; j++){
    xgrid[j] = mcurve.quad[0][j] * 0.5 * mcurve.size;
    ygrid[j] = mcurve.quad[1][j] * 0.5 * mcurve.size;
  }

  make_morton_grid(mcurve, xgrid, ygrid, 4);
  make_morton_key(&mcurve);

/* ---------------- Write on File ------------------- */
  out = fopen("morton_curve.dat", "w");

  for(i = 0; i < mcurve.npoints; i++)
    //fprintf(out, "%lg\t%lg\n", mcurve.xpos[i], mcurve.ypos[i]);
    fprintf(out, "%lg\t%lg\t%d\n", mcurve.xpos[i], mcurve.ypos[i], mcurve.key_list[i]);

  fclose(out);
/* ------------------------------------------------- */
// Deallocate
  morton_dealloc(mcurve);
  delete [] xgrid;
  delete [] ygrid;

  return 0;
}
