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

typedef std::bitset<64> bittype; // Alias

bool mybitf(bittype & n, int pos){
  return (n[pos] == 1);
}

void mybitset(bittype & n, int pos){
  n[pos] = 1;
}

void morton_lattice_init(struct Morton * mort, int scan, int size){

  int tot_num_of_points = (int) std::pow(4, scan);
  double * tmp = new double[8]; // This is equivalent to malloc(8 * sizeof(double)) in C

  (* mort).npoints = tot_num_of_points;
  //mort->npoints = tot_num_of_points;
  //mort.npoints = tot_num_of_points; se metto &
  (* mort).xpos = new double [tot_num_of_points];
  (* mort).ypos = new double [tot_num_of_points];
  (* mort).key_list = new int [tot_num_of_points];


  (* mort).reclevel = 0; // recursion level
  (* mort).recdepth = scan; // recursion depth
  (* mort).size = size;

  // Array of unit box coordinates
  (* mort).quad = new double * [2];
  tmp = new double[8];

  (* mort).quad[0] = &tmp[0]; // here we will store x values
  (* mort).quad[1] = &tmp[4]; // here we will store y values

  (* mort).quad[0][0] = 0.5;
  (* mort).quad[1][0] = 0.5;

  (* mort).quad[0][1] = 1.5;
  (* mort).quad[1][1] = 0.5;

  (* mort).quad[0][2] = 0.5;
  (* mort).quad[1][2] = 1.5;

  (* mort).quad[0][3] = 1.5;
  (* mort).quad[1][3] = 1.5;

  // Array of coordinates of the unit box corners
  (* mort).corner = new double * [2];
  tmp = new double[8];

  (* mort).corner[0] = &tmp[0]; // here we will store x values
  (* mort).corner[1] = &tmp[4]; // here we will store y values

  (* mort).corner[0][0] = 0.0;
  (* mort).corner[1][0] = 0.0;

  (* mort).corner[0][1] = 1.0;
  (* mort).corner[1][1] = 0.0;

  (* mort).corner[0][2] = 0.0;
  (* mort).corner[1][2] = 1.0;

  (* mort).corner[0][3] = 1.0;
  (* mort).corner[1][3] = 1.0;

}

void morton_dealloc(struct Morton * mort){

  delete [] (* mort).xpos;
  delete [] (* mort).ypos;
  delete [] (* mort).key_list;
  delete [] (* mort).quad[0];
  delete [] (* mort).corner[0];
  delete [] (* mort).quad;
  delete [] (* mort).corner;
}

void make_morton_grid(struct Morton * mort,
                      double * xgrid,
                      double * ygrid,
                      int n){

  int i, j, iad;

  // reclevel starts from 0
  if( (* mort).reclevel + 1 >= (* mort).recdepth ){
    /* COPY RESULTS INTO X_GRID AND Y_GRID */

    printf("\n\tNpoints = %d, n = %d", (* mort).npoints, n);

    for(i = 0; i < n; i++){
      (* mort).xpos[i] = xgrid[i];
      (* mort).ypos[i] = ygrid[i];
    }

    return;
  }

  double* xsub = new double[4 * n];
  double* ysub = new double[4 * n];

  (*mort).reclevel += 1;

  for(i = 0; i < 4; i++){

    iad = i * n;

    for(j = 0; j < n; j++){
      xsub[iad + j] = (xgrid[j] + (* mort).corner[0][i] * (* mort).size ) * 0.5;
      ysub[iad + j] = (ygrid[j] + (* mort).corner[1][i] * (* mort).size ) * 0.5;
    }
  }
  make_morton_grid(mort, xsub, ysub, 4*n);
  delete [] xsub;
  delete [] ysub;
}


void make_morton_key(struct Morton* mort){

  int ix, iy, tmp;
  bittype ixb, iyb, tmpb;
  int j, levmorton = 8, levkey = 0;

  double cj = std::pow(2, levmorton) / (* mort).size;


  for(j = 0; j < (* mort).npoints; j++){

    ix = (* mort).xpos[j] * cj;
    //ix = (* mort).xpos[j];
    iy = (* mort).ypos[j] * cj;
    //iy = (* mort).ypos[j];
    tmp = 0;

    ixb = bittype(ix);
    iyb = bittype(iy);
    tmpb = bittype(tmp);

    while(levkey < levmorton){

      if(mybitf(ixb, levkey))
	mybitset(tmpb, 2*levkey);

      if(mybitf(iyb, levkey))
	mybitset(tmpb, 2*levkey+1);

      levkey++;
    }
    (* mort).key_list[j] = tmpb.to_ulong();
    levkey = 0;
  }

}


int main(){

  int i, j;
  double * xgrid, * ygrid;
  FILE* out;
  struct Morton mcurve;

  xgrid = new double [4];
  ygrid = new double [4];

  morton_lattice_init(&mcurve, 4, 16);

  printf("\n\tmorton curve lscan: %d\n\n", mcurve.recdepth);
  for(i = 0; i < 2; i++){
    for(j = 0; j < 4; j++)
    printf("\t%lg", mcurve.quad[i][j]);
    printf("\n");
  }

  /* initializing xgrid and ygrid */
  for(j = 0; j < 4; j++){
    xgrid[j] = mcurve.quad[0][j] * 0.5 * mcurve.size;
    ygrid[j] = mcurve.quad[1][j] * 0.5 * mcurve.size;
  }

  make_morton_grid(&mcurve, xgrid, ygrid, 4);
  make_morton_key(&mcurve);

  //std::cout << "\n\n\tsaving output\n" << std::endl;
  printf("\n\n\tsaving output\n");
  out = fopen("morton_curve.dat", "w");

  for(i = 0; i < mcurve.npoints; i++)
    //fprintf(out, "%lg\t%lg\n", mcurve.xpos[i], mcurve.ypos[i]);
    fprintf(out, "%lg\t%lg\t%d\n", mcurve.xpos[i], mcurve.ypos[i], mcurve.key_list[i]);

  fclose(out);
  morton_dealloc(&mcurve);

  delete [] xgrid;
  delete [] ygrid;

  return 0;
}
