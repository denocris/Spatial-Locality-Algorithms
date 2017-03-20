#include "tools.h"


int main(){

  int i, j, niterations;
  double * xgrid, * ygrid;
  FILE * out;
  struct Hilbert hcurve;

  xgrid = new double [4];
  ygrid = new double [4];

  niterations = 3;

  hilbert_lattice_init(hcurve, niterations, 16);

/* ------------------- Initialization ---------------- */
  for(j = 0; j < 4; j++){
    xgrid[j] = hcurve.quad_pos[0][j] * 0.5 * hcurve.side;
    ygrid[j] = hcurve.quad_pos[1][j] * 0.5 * hcurve.side;
  }

  make_hilbert_grid(hcurve, xgrid, ygrid, 4);

/* ---------------- Write on File ------------------- */
  out = fopen("hilbert_curve.dat", "w");

  for(i = 0; i < hcurve.npoints; i++)
    fprintf(out, "%lg\t%lg\n", hcurve.xpos[i], hcurve.ypos[i]);

  fclose(out);
/* ------------------------------------------------- */
// Deallocate
  hilbert_dealloc(hcurve);
  delete [] xgrid;
  delete [] ygrid;

  return 0;
}
