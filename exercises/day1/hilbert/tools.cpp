#include "tools.h"

void hilbert_lattice_init(struct Hilbert & hilb, int scan, int side){

  int i, j, k;
  double * tmp = new double [8];

  int tot_num_of_points;
  tot_num_of_points = std::pow(4, scan);

  hilb.npoints = tot_num_of_points;

  hilb.side = side;
  hilb.recdepth = scan;
  hilb.reclevel = 0;

  hilb.xpos = new double [tot_num_of_points];
  hilb.ypos = new double [tot_num_of_points];


  hilb.quad_pos = new double * [2];
  hilb.corner = new double * [2];

  hilb.quad_pos[0] = &tmp[0];
  hilb.quad_pos[1] = &tmp[4];

  hilb.quad_pos[0][0] = 0.5;
  hilb.quad_pos[1][0] = 0.5;

  hilb.quad_pos[0][1] = 0.5;
  hilb.quad_pos[1][1] = 1.5;

  hilb.quad_pos[0][2] = 1.5;
  hilb.quad_pos[1][2] = 1.5;

  hilb.quad_pos[0][3] = 1.5;
  hilb.quad_pos[1][3] = 0.5;

  tmp = new double [8];

  hilb.corner[0] = &tmp[0];
  hilb.corner[1] = &tmp[4];

  hilb.corner[0][0] = 0.0;
  hilb.corner[1][0] = 0.0;

  hilb.corner[0][1] = 0.0;
  hilb.corner[1][1] = 1.0;

  hilb.corner[0][2] = 1.0;
  hilb.corner[1][2] = 1.0;

  hilb.corner[0][3] = 1.0;
  hilb.corner[1][3] = 0.0;

}


void hilbert_dealloc(struct Hilbert & hilb){

  delete [] hilb.xpos;
  delete [] hilb.ypos;
  delete [] hilb.corner;
  delete [] hilb.quad_pos;
}


void make_hilbert_grid(struct Hilbert & hilb,
                      double * xgrid,
                      double * ygrid,
                      int n){

  int i, j, iad, rot;

  if( hilb.reclevel + 1 >= hilb.recdepth ){

    printf("\n\tNpoints = %d, n = %d\n", hilb.npoints, n);

    for(i = 0; i < n; i++){
      hilb.xpos[i] = xgrid[i];
      hilb.ypos[i] = ygrid[i];
    }

    return;
  }

  double * xsub = new double[4 * n];
  double * ysub = new double[4 * n];

  double * xsubrr = new double[4 * n];
  double * ysubrr = new double[4 * n];

  hilb.reclevel += 1;

  for(i = 0; i < 4; i++){ // Copy in each subq

    if(i == 0){
      iad = i * n;

      for(j = 0; j < n; j++){
        xsub[iad + j] = (ygrid[j] + hilb.corner[0][i] * hilb.side ) * 0.5;
        ysub[iad + j] = (xgrid[n - 1 - j] + hilb.corner[1][i] * hilb.side ) * 0.5;

        }

      for(j = 0; j < n; j++){
        xsubrr[iad + j] = xsub[iad + j];
        ysubrr[iad + j] = ysub[iad + (n - 1 - j)];
        }
    }

    else if(i == 3){
      iad = i * n;

      for(j = 0; j < n; j++){
        int tmp = hilb.side - ygrid[n - 1 - j];
        xsub[iad + j] = ( tmp + hilb.corner[0][i] * hilb.side ) * 0.5;
        ysub[iad + j] = (xgrid[j] + hilb.corner[1][i] * hilb.side ) * 0.5;
        }

      for(j = 0; j < n; j++){
        xsubrr[iad + j] = xsub[iad + j];
        ysubrr[iad + j] = ysub[iad + (n - 1 - j)];
        }
    }

    else if(i == 1 || i ==2){
      iad = i * n;

      for(j = 0; j < n; j++){
        xsubrr[iad + j] = (xgrid[j] + hilb.corner[0][i] * hilb.side ) * 0.5;
        ysubrr[iad + j] = (ygrid[j] + hilb.corner[1][i] * hilb.side ) * 0.5;
        }
    }


  }

  make_hilbert_grid(hilb, xsubrr, ysubrr, 4*n);
  delete [] xsub;
  delete [] xsubrr;
  delete [] ysub;
  delete [] ysubrr;
}
