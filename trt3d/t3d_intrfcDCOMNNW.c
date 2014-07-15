extern int NNWSETELMENT(void);
extern void calculate_d1d2d3(int, int, double, double, double, double, double, double[], double[], double[], double[], double[4]);

void intrfcdcomnnw_(int *seed_ptr, int *num_ions_ptr, double *rho_ptr, double *Er_ptr, 
	   double *RAX_ptr, double *beta_ptr, double *BB_ptr, 
	   double TT[3], double NN[3], double ZZ[3], double MM[3], 
	   double DD[4])
{
  int seed;
  int num_ions;
  double rho;
  double Er;
  double RAX;
  double beta;
  double BB;


  NNWSETELMENT();	 /*  Declaration of using the NNW data base.*/

  seed = *seed_ptr;
  num_ions = *num_ions_ptr;
  Er = *Er_ptr;
  rho = *rho_ptr;
  RAX = *RAX_ptr;
  beta = *beta_ptr;
  BB = *BB_ptr;
//    if(rho<0.1){
//        calculate_d1d2d3(seed, num_ions, 0.1, 0.0, RAX, beta, BB, TT, NN, ZZ, MM, DD);
//    }
//    else {
        calculate_d1d2d3(seed, num_ions, rho, Er, RAX, beta, BB, TT, NN, ZZ, MM, DD);
//    }

//  printf("%d,%d\n%e,%e,%e,%e,%e\n%e,%e,%e,\n%e,%e,%e\n%e,%e,%e\n%e,%e,%e\n%e,%e,%e,%e\n",seed, num_ions, rho, Er, RAX, beta, BB, TT[0],TT[1],TT[2], NN[0],NN[1],NN[2], ZZ[0],ZZ[1],ZZ[2], MM[0],MM[1],MM[2], DD[0],DD[1],DD[2],DD[3]);
//  calculate_d1d2d3(seed, num_ions, rho, Er, RAX, beta, BB, TT, NN, ZZ, MM, DD);
}
