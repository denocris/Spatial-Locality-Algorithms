#include <string>
#include <iostream>
#include <cmath>
#include <bitset>


typedef std::bitset<64> bittype; // Alias

typedef struct Hilbert{
  int npoints;
  double * xpos;
  double * ypos;
  double **  quad_pos;
  double **  corner;
  int reclevel, recdepth, side;

} hcurve;



void hilbert_lattice_init(struct Hilbert * hilb, int scan, int side){

  int i, j, k;
  double * tmp = new double [8];
  //double** tr = new double * [8];

  int tot_num_of_points;
  tot_num_of_points = std::pow(4, scan);

  (* hilb).npoints = tot_num_of_points;

  (* hilb).side = side;
  (* hilb).recdepth = scan;
  (* hilb).reclevel = 0;

  (* hilb).xpos = new double [tot_num_of_points];
  (* hilb).ypos = new double [tot_num_of_points];



  (* hilb).quad_pos = new double * [2];
  (* hilb).corner = new double * [2];

  (* hilb).quad_pos[0] = &tmp[0];
  (* hilb).quad_pos[1] = &tmp[4];

  (* hilb).quad_pos[0][0] = 0.5;
  (* hilb).quad_pos[1][0] = 0.5;

  (* hilb).quad_pos[0][1] = 0.5;
  (* hilb).quad_pos[1][1] = 1.5;

  (* hilb).quad_pos[0][2] = 1.5;
  (* hilb).quad_pos[1][2] = 1.5;

  (* hilb).quad_pos[0][3] = 1.5;
  (* hilb).quad_pos[1][3] = 0.5;

  tmp = new double [8];

  (* hilb).corner[0] = &tmp[0];
  (* hilb).corner[1] = &tmp[4];

  (* hilb).corner[0][0] = 0.0;
  (* hilb).corner[1][0] = 0.0;

  (* hilb).corner[0][1] = 0.0;
  (* hilb).corner[1][1] = 1.0;

  (* hilb).corner[0][2] = 1.0;
  (* hilb).corner[1][2] = 1.0;

  (* hilb).corner[0][3] = 1.0;
  (* hilb).corner[1][3] = 0.0;



}



void hilbert_dealloc(struct Hilbert * hilb){

  delete [] (* hilb).xpos;
  delete [] (* hilb).ypos;
  delete [] (* hilb).corner;
  delete [] (* hilb).quad_pos;
}


void make_hilbert_grid(struct Hilbert * hilb,
                      double * xgrid,
                      double * ygrid,
                      int n){

  int i, j, iad, rot;

  if( (* hilb).reclevel + 1 >= (* hilb).recdepth ){

    printf("\n\tNpoints = %d, n = %d", (* hilb).npoints, n);

    for(i = 0; i < n; i++){
      (* hilb).xpos[i] = xgrid[i];
      (* hilb).ypos[i] = ygrid[i];
    }

    return;
  }

  double * xsub = new double[4 * n];
  double * ysub = new double[4 * n];

  double * xsubr = new double[4 * n];
  double * ysubr = new double[4 * n];

  double * xsubrr = new double[4 * n];
  double * ysubrr = new double[4 * n];

  (* hilb).reclevel += 1;

  for(i = 0; i < 4; i++){ // Copy in each subq

    iad = i * n;

    for(j = 0; j < n; j++){
      xsub[iad + j] = (xgrid[j] + (* hilb).corner[0][i] * (* hilb).side ) * 0.5;
      ysub[iad + j] = (ygrid[j] + (* hilb).corner[1][i] * (* hilb).side ) * 0.5;
    }
  }

  
  for(i = 0; i < 4; i++){ // Rotation and reflection in subq 0 and 3


    if(i == 0 || i == 3){
      iad = i * n;

      if(i == 0){rot = 1;}
      if(i == 3){rot = - 1;}


      for(j = 0; j < n; j++){  // Rotation subq 0
          xsubr[iad + j] = xsub[iad + (j + rot) % n];
          ysubr[iad + j] = ysub[iad + (((j + rot) % n)+n)%n];
        }

      for(j = 0; j < n; j++){ // Reflection subq 0
        xsubrr[iad + j] = xsubr[iad + (n - 1 - j)];
        ysubrr[iad + j] = ysubr[iad + (n - 1 - j)];
      }
  }

    else if(i == 1 || i == 2) {
      iad = i * n;
        for(j = 0; j < n; j++){
            xsubr[iad + j] = xsub[iad + j];
            ysubr[iad + j] = ysub[iad + j];
            xsubrr[iad + j] = xsub[iad + j];
            ysubrr[iad + j] = ysub[iad + j];
          }

  }

  }

  // for(j = 0; j < 4*n; j++) printf("%lg", xsub[j] );
  // printf("\n");
  // for(j = 0; j < 4*n; j++) printf("%lg", xsubr[j] );
  // printf("\n");
  // for(j = 0; j < 4*n; j++) printf("%lg", xsubrr[j] );
  // printf("\n");
  // printf("-------\n");
  // for(j = 0; j < 4*n; j++) printf("%lg", ysub[j] );
  // printf("\n");
  // for(j = 0; j < 4*n; j++) printf("%lg", ysubr[j] );
  // printf("\n");
  // for(j = 0; j < 4*n; j++) printf("%lg", ysubrr[j] );



  make_hilbert_grid(hilb, xsubrr, ysubrr, 4*n);
  //make_hilbert_grid(hilb, xsub, ysub, 4*n);
  delete [] xsub;
  delete [] xsubr;
  delete [] xsubrr;
  delete [] ysub;
  delete [] ysubr;
  delete [] ysubrr;
}



int main(){

  int i, j;
  double * xgrid, * ygrid;
  FILE* out;
  struct Hilbert hcurve;

  xgrid = new double [4];
  ygrid = new double [4];

  hilbert_lattice_init(&hcurve, 2, 16);

  printf("\n\thilbert curve lscan: %d\n\n", hcurve.recdepth);
  for(i = 0; i < 2; i++){
    for(j = 0; j < 4; j++)
    printf("\t%lg", hcurve.quad_pos[i][j]);
    printf("\n");
  }

  for(j = 0; j < 4; j++){
    xgrid[j] = hcurve.quad_pos[0][j] * 0.5 * hcurve.side;
    ygrid[j] = hcurve.quad_pos[1][j] * 0.5 * hcurve.side;
    //printf("x %lg\n", xgrid[j]);
    //printf("y %lg\n", ygrid[j]);
  }

  make_hilbert_grid(&hcurve, xgrid, ygrid, 4);

  //std::cout << "\n\n\tsaving output\n" << std::endl;
  printf("\n\n\tsaving output\n");
  out = fopen("hilbert_curve.dat", "w");

  for(i = 0; i < hcurve.npoints; i++)
    fprintf(out, "%lg\t%lg\n", hcurve.xpos[i], hcurve.ypos[i]);

  fclose(out);
  hilbert_dealloc(&hcurve);
  delete [] xgrid;
  delete [] ygrid;

  return 0;
}
