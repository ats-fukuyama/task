/*  NNW_DATA.c	version. 20071012 */
/*  更新記録 2007.10.12 有限ベータをNNWに組み込み（WV[15]シリーズ） */
#define nu_star_extrapolation 0	/* .TURE. 1 or .FALSE. 0 || nu_starの外装を行うかどうか */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "nandTprofile.DCOM_NNW.h"

#define NUM_IONS_TMP 10
/*// Declaration of fundamental parameter about NNW //*/
#ifndef use_NNWwithBETA
#define NUM_INPUT 3
#endif
#ifdef use_NNWwithBETA
#define NUM_INPUT 4
#endif
#define NUM_HIDDEN  15
#define NUM_OUTPUT 1
#define NUM_OMEGA ((NUM_INPUT+1)*NUM_HIDDEN+(NUM_HIDDEN+1)*NUM_OUTPUT)
#define NUM_RAX 20

//#define OLD_NNW_20091009
#define NEW_NNW_20121019


/*// Declaration of structures //*/
typedef struct Weightvalue{
	double w1[NUM_INPUT+1][NUM_HIDDEN+1];
	double w2[NUM_HIDDEN+1][NUM_OUTPUT];
	double omega[NUM_OMEGA];
} WEIGHTVALUE;
typedef struct NNWvalue{
	double xunit[NUM_INPUT+1];
	double hiddenunit[NUM_HIDDEN+1];
	double yunit[NUM_OUTPUT+1];
} NNWVALUE;

/*// Declaration of global variables //*/
WEIGHTVALUE WV[NUM_RAX];
double nu_star_MAX[NUM_RAX];
double nu_star_MIN[NUM_RAX];
double D_star_MAX[NUM_RAX];
double D_star_MIN[NUM_RAX];
double GGGMAX[NUM_RAX];

/*// Declaration of functions //*/
int NNWSETELMENT(void);
double CALCULATE_DIFFUSION_COEFFICIENT(double, double, double, double, double);
int ScalingInput(double[NUM_INPUT], double, double, double, double, double);
int CALCULATEbyNNW(struct Weightvalue*, struct NNWvalue*, double, double);
double DERF(double);
double CLOG(double, double);
double ERRORFUNCTION(double);
void calculate_d1d2d3(int, int, double, double, double, double, double, double[], double[], double[], double[], double[4]);
void SET_ARADI_IOTA(double, double, double, double*, double*);

void calculate_d1d2d3(int seed, int num_ions, double rho, double Er, double RAX, double beta,
					 double BB, double TT[], double NN[], double ZZ[], double MM[],
					 double DD[])
{

	#define  vth_times 3
	#define  dvpara 100
	#define MAX_iv (vth_times * dvpara)
	double CHARGE =(1.60217733e-19);
	double CC=  (2.99792458e+8);
	double PI=	(3.141592653589793e0);
	double MASSE=	(9.1093897e-31);
	double MASSI=	(1.67262305e-27);
	double epsilon0=(1.0e7/4.0/PI/CC/CC);

	int i;
	int iv, itarget;
	double vth[NUM_IONS_TMP+1]={0.0}, herf[NUM_IONS_TMP+1]={0.0}, clog[NUM_IONS_TMP+1]={0.0};
	double v, dv, GGG, D_star, Dk, D_plateau[NUM_IONS_TMP+1]={0.0};
	double nu_origin, nu_origin_tmp[NUM_IONS_TMP+1]={0.0};
	double D1_tmp[MAX_iv+1]={0.0}, D2_tmp[MAX_iv+1]={0.0}, D3_tmp[MAX_iv+1]={0.0};
	double D1_SUM =0.0,  D2_SUM =0.0,  D3_SUM =0.0;
	double nu_star;
	double mass[NUM_IONS_TMP+1];
	double aradi, iota;

	/* setteing of aradi and iota */
	SET_ARADI_IOTA(RAX, beta, rho, &aradi, &iota);

	/* setteing of the mass */
	for(i=0;i<num_ions+1;i++){
		if(MM[i]<0){
			mass[i]=MASSE;
		}
		else{
			mass[i]=MM[i]*MASSI;
		}
	}

	for(i=0;i<num_ions+1;i++){
		vth[i] = sqrt((double)2.0*TT[i]*CHARGE/mass[i]);
		D_plateau[i] = (PI)*pow(vth[i],3)*pow(mass[i],2)/iota
						/pow((ZZ[i]*CHARGE),2)/pow(BB,2)/RAX/(double)16.0;
//        printf("%e,%e,%e,%e\n",vth[i],mass[i],iota,ZZ[i]);        
	}
	
/* vth 表示 */
/*
		for(itarget=0;itarget<num_ions+1;itarget++){
			herf[itarget]=ERRORFUNCTION((vth[seed]/vth[itarget]));

			if(NN[itarget]<1e-300){
				clog[itarget]=0.0;
			}
			else
			{
			   clog[itarget] = (double)24.0-log(sqrt(NN[0]*(double)1.e-6)/TT[0]);
			}
		}

		for(itarget=0;itarget<num_ions+1;itarget++){
			nu_origin_tmp[itarget] = pow((ZZ[itarget]*CHARGE),2)*pow((ZZ[seed]*CHARGE),2)*NN[itarget]*clog[itarget]
								/(double)4.0/PI/pow(epsilon0,2)/pow(mass[seed],2)*herf[itarget]/(vth[seed]*vth[seed]*vth[seed]);
								//printf("seed: %d nu_origin_tmp[%d]: %e clog: %e herf: %e\n", seed, itarget, nu_origin_tmp[itarget], clog[itarget], herf[itarget]);
		}

		nu_origin =0.0;
		for(itarget=0;itarget<num_ions+1;itarget++){
			nu_origin =nu_origin + nu_origin_tmp[itarget];
		}
		nu_star = RAX*nu_origin/iota/vth[seed];
		GGG = RAX*Er/iota/(rho*aradi)/vth[seed]/BB;
		D_star = CALCULATE_DIFFUSION_COEFFICIENT(nu_star, rho, GGG, RAX, beta);
//		Dk =  D_plateau[seed]*D_star *(pow((v/vth[seed]),3));
		printf("Rax: %.2f TT: %e seed: %d nu_oringin:%e nu_star: %e GGG: %e Dplateau :%e Dstar: %e \n", RAX, TT[seed]/1e3, seed, nu_origin, nu_star, GGG, D_plateau[seed], D_star);
        exit(0);
*/
		/* ここまで */


	dv = vth[seed] /(double)dvpara;

	for (iv=1;iv<MAX_iv+1;iv++){ /*/ **** dv loop **** /*/

		v=(double)iv*dv;

		GGG = RAX*Er/iota/(rho*aradi)/v/BB;
		for(itarget=0;itarget<num_ions+1;itarget++){
			herf[itarget]=ERRORFUNCTION((v/vth[itarget]));
			if(NN[itarget]<1e-300){
				clog[itarget]=0.0;
			}
			else
			{
			   clog[itarget] = (double)24.0-log(sqrt(NN[0]*(double)1.e-6)/TT[0]);
			   /* OLD clog */	
#ifdef OLDNUandCLOG
			  clog[itarget] = log((double)12.0*PI*pow(sqrt(epsilon0),3)/CHARGE/CHARGE/CHARGE) + log(pow(sqrt(TT[itarget]*CHARGE),3)/sqrt(NN[itarget]));
#endif
			}
		}

		for(itarget=0;itarget<num_ions+1;itarget++){
			nu_origin_tmp[itarget] = pow((ZZ[itarget]*CHARGE),2)*pow((ZZ[seed]*CHARGE),2)*NN[itarget]*clog[itarget]
								/(double)4.0/PI/pow(epsilon0,2)/pow(mass[seed],2)*herf[itarget]/(v*v*v);
			/* OLD collision frequency */
#ifdef OLDNUandCLOG
			nu_origin_tmp[itarget] = pow((ZZ[itarget]*CHARGE),2)*pow((ZZ[seed]*CHARGE),2)*NN[itarget]*clog[itarget]
							/(double)2.0/PI/pow(epsilon0,2)/pow(mass[seed],2)*herf[itarget]/(v*v*v);
#endif
		}

		nu_origin =0.0;
		for(itarget=0;itarget<num_ions+1;itarget++){
			nu_origin =nu_origin + nu_origin_tmp[itarget];
		}
		nu_star = RAX*nu_origin/iota/v;
		//printf("vth[0] = %e\tvth[1] = %e\tv %e:iota %e:nu_star %e\n", vth[0], vth[1],v,iota, nu_star);/* V_TERMAL_TEST!!! 20071015  */
		//system("pause");
//                if(betaip_flg==0){

        D_star = CALCULATE_DIFFUSION_COEFFICIENT(nu_star, rho, GGG, RAX, beta);
        
		Dk =  D_plateau[seed]*D_star *(pow((v/vth[seed]),3));
		D1_tmp[iv] =  (Dk  * (pow((v /vth[seed]),2))*exp((double)-1*(pow((v /vth[seed]),2)))/(vth[seed]));
		D2_tmp[iv] =  (Dk  * (pow((v /vth[seed]),4))*exp((double)-1*(pow((v /vth[seed]),2)))/(vth[seed]));
		D3_tmp[iv] =  (Dk  * (pow((v /vth[seed]),6))*exp((double)-1*(pow((v /vth[seed]),2)))/(vth[seed]));
//                }
//              else if(betaip_flg==1){
//                double D_star_low,D_star_high;
//		D_star_low = CALCULATE_DIFFUSION_COEFFICIENT(nu_star, rho, GGG, RAX, beta_low);
//		D_star_high = CALCULATE_DIFFUSION_COEFFICIENT(nu_star, rho, GGG, RAX, beta_high);
//               D_star=(D_star_high-D_star_low)/(beta_high-beta_low)*beta +D_star_low;
//		Dk =  D_plateau[seed]*D_star *(pow((v/vth[seed]),3));
//		D1_tmp[iv] =  (Dk  * (pow((v /vth[seed]),2))*exp((double)-1*(pow((v /vth[seed]),2)))/(vth[seed]));
//		D2_tmp[iv] =  (Dk  * (pow((v /vth[seed]),4))*exp((double)-1*(pow((v /vth[seed]),2)))/(vth[seed]));
//		D3_tmp[iv] =  (Dk  * (pow((v /vth[seed]),6))*exp((double)-1*(pow((v /vth[seed]),2)))/(vth[seed]));
//                }
	}
	for (iv=1;iv<MAX_iv-3;iv++){
		/* Newton-Cotes formulas */
		D1_SUM = D1_SUM + (D1_tmp[iv]+ 3.0*D1_tmp[iv+1]+ 3.0*D1_tmp[iv+2]+ D1_tmp[iv+3])/(double)8.0;
//		printf("%d\t%e\n", iv, D1_SUM );
		D2_SUM = D2_SUM + (D2_tmp[iv]+ 3.0*D2_tmp[iv+1]+ 3.0*D2_tmp[iv+2]+ D2_tmp[iv+3])/(double)8.0;
		D3_SUM = D3_SUM + (D3_tmp[iv]+ 3.0*D3_tmp[iv+1]+ 3.0*D3_tmp[iv+2]+ D3_tmp[iv+3])/(double)8.0;		
	} 
//	exit(0);
	DD[1]=(double)4/sqrt((double)PI)*D1_SUM*dv;
	DD[2]=(double)4/sqrt((double)PI)*D2_SUM*dv;
	DD[3]=(double)4/sqrt((double)PI)*D3_SUM*dv;

	return;
}

double CALCULATE_DIFFUSION_COEFFICIENT(double nu_star, double rho, double GGG, double RAX, double beta)
{

	double D_star;

	NNWVALUE NNWV;
	int iRAX;

	iRAX=(int)(RAX * 100.0);

	switch(iRAX){
	case 345:if(beta==0.0){iRAX = 0;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 350:if(beta==0.0){iRAX = 1;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 353:if(beta==0.0){iRAX = 2;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 360:if(beta==0.0){iRAX = 3;break;}else 
        if(beta!=0.0){printf("BAD CONDITION\n");exit(1);}else
        if(fabs(beta-0.01)<1e-5) {iRAX = 6;break;}else 
            if(fabs(beta-0.005)<1e-5){iRAX = 9;break;}else 
                if(fabs(beta-0.02)<1e-5) {iRAX = 7;break;}else 
                    if(beta==0.03) {iRAX = 8;break;}else 
			 {printf("(1)An improper value was substituted for beta. \n");exit(1);}
	case 375:
#ifndef use_NNWwithBETA
			 if(beta==0.0){iRAX = 4;break;}else
#endif
#ifdef use_NNWwithBETA
			 if(beta>=0.00 && beta<=0.05){iRAX = 15;break;}else
#endif
			 {printf("An improper value was substituted for beta. \n");exit(1);}
	case 390:if(beta==0.0){iRAX = 5;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	default :printf("An improper value was substituted for RAX. \n");exit(1);
	}

	if (nu_star==0){
		D_star=0.0;
		return(D_star);
	}
	else if (nu_star <(nu_star_MIN[iRAX])){
#if (nu_star_extrapolation==0)
		nu_star = (double)nu_star_MIN[iRAX];
#endif
#if (nu_star_extrapolation!=0)
		double GGG0, rho0, RAX0, beta0;
		double nu_star_1, nu_star_2;
		double  D_star_1, D_star_2;

		if(GGG<0){
			if(GGG<(-1.0*GGGMAX[iRAX])){
				GGG=-1.0*GGGMAX[iRAX];
			}
			GGG = GGG * -(double)1.0;
		}
		else{
			if(GGG>GGGMAX[iRAX]){
				GGG=GGGMAX[iRAX];
			}
		}

		GGG0=GGG;
		rho0=rho;
		RAX0=RAX;
		beta0=beta;
		nu_star_1=(double)nu_star_MIN[iRAX];

		ScalingInput(NNWV.xunit, nu_star_1, GGG, rho, RAX, beta);
		CALCULATEbyNNW(&WV[iRAX], &NNWV, RAX, beta);
		D_star_1 = exp((((NNWV.yunit[0]+(double)0.9)/(double)1.8)*(log(D_star_MAX[iRAX])-log(D_star_MIN[iRAX])))+log(D_star_MIN[iRAX]));
	
		GGG=GGG0;
		rho=rho0;
		RAX=RAX0;
		beta=beta0;
		nu_star_2=(double)nu_star_MIN[iRAX]*10.0;

		ScalingInput(NNWV.xunit, nu_star_2, GGG, rho, RAX, beta);
		CALCULATEbyNNW(&WV[iRAX], &NNWV, RAX, beta);
		D_star_2 = exp((((NNWV.yunit[0]+(double)0.9)/(double)1.8)*(log(D_star_MAX[iRAX])-log(D_star_MIN[iRAX])))+log(D_star_MIN[iRAX]));

		D_star=((log10(D_star_1)-log10(D_star_2))/(log10(nu_star_1)-log10(nu_star_2)))*(log10(nu_star)-log10(nu_star_1))+log10(D_star_1);
		D_star=pow(10, D_star);

		return(D_star);
#endif
	}
	else if(nu_star > nu_star_MAX[iRAX]){
#if(nu_star_extrapolation==0)
		nu_star = (double)nu_star_MAX[iRAX];
#endif
#if (nu_star_extrapolation!=0)
		double GGG0, rho0, RAX0, beta0;
		double nu_star_1, nu_star_2;
		double  D_star_1, D_star_2;

		if(GGG<0){
			if(GGG<(-1.0*GGGMAX[iRAX])){
				GGG=-1.0*GGGMAX[iRAX];
			}
			GGG = GGG * -(double)1.0;
		}
		else{
			if(GGG>GGGMAX[iRAX]){
				GGG=GGGMAX[iRAX];
			}
		}

		GGG0=GGG;
		rho0=rho;
		RAX0=RAX;
		beta0=beta;
		nu_star_1=(double)nu_star_MAX[iRAX];

		ScalingInput(NNWV.xunit, nu_star_1, GGG, rho, RAX, beta);
		CALCULATEbyNNW(&WV[iRAX], &NNWV, RAX, beta);
		D_star_1 = exp((((NNWV.yunit[0]+(double)0.9)/(double)1.8)*(log(D_star_MAX[iRAX])-log(D_star_MIN[iRAX])))+log(D_star_MIN[iRAX]));
	
		GGG=GGG0;
		rho=rho0;
		RAX=RAX0;
		beta=beta0;
		nu_star_2=(double)nu_star_MAX[iRAX]/10.0;

		ScalingInput(NNWV.xunit, nu_star_2, GGG, rho, RAX, beta);
		CALCULATEbyNNW(&WV[iRAX], &NNWV, RAX, beta);
		D_star_2 = exp((((NNWV.yunit[0]+(double)0.9)/(double)1.8)*(log(D_star_MAX[iRAX])-log(D_star_MIN[iRAX])))+log(D_star_MIN[iRAX]));

		D_star=((log10(D_star_1)-log10(D_star_2))/(log10(nu_star_1)-log10(nu_star_2)))*(log10(nu_star)-log10(nu_star_1))+log10(D_star_1);
		D_star=pow(10, D_star);

		return(D_star);
#endif
	}



	if(GGG<0){
		if(GGG<(-1.0*GGGMAX[iRAX])){
			GGG=-1.0*GGGMAX[iRAX];
		}
		GGG = GGG * -(double)1.0;
	}
	else{
		if(GGG>GGGMAX[iRAX]){
			GGG=GGGMAX[iRAX];
		}
	}
		
	ScalingInput(NNWV.xunit, nu_star, GGG, rho, RAX, beta);

	CALCULATEbyNNW(&WV[iRAX], &NNWV, RAX, beta);
	
	D_star = exp((((NNWV.yunit[0]+(double)0.9)/(double)1.8)*(log(D_star_MAX[iRAX])-log(D_star_MIN[iRAX])))+log(D_star_MIN[iRAX]));

	return(D_star);
}

int ScalingInput(double x[NUM_INPUT], double nu_star, double GGG, double rho, double RAX, double beta){
	int iRAX;

	iRAX=(int)(RAX * 100);

	switch(iRAX){
	case 345:if(beta==0.0){iRAX = 0;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 350:if(beta==0.0){iRAX = 1;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 353:if(beta==0.0){iRAX = 2;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 360:if(beta==0.0){iRAX = 3;break;}else 
        if(fabs(beta-0.01)<1e-5) {iRAX = 6;break;}else 
            if(fabs(beta-0.005)<1e-5){iRAX = 9;break;}else 
                if(fabs(beta-0.02)<1e-5) {iRAX = 7;break;}else 
                    if(beta==0.03) {iRAX = 8;break;}else 
			 {printf("(2)An improper value was substituted for beta. \n");exit(1);}
	case 375:
#ifndef use_NNWwithBETA
			 if(beta==0.0){iRAX = 4;break;}else
#endif
#ifdef use_NNWwithBETA
			 if(beta>=0.00 && beta<=0.05){iRAX = 15;break;}else
#endif
			 {printf("An improper value was substituted for beta. \n");exit(1);}
	case 390:if(beta==0.0){iRAX = 5;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	default :printf("An improper value was substituted for RAX. \n");exit(1);
	}


	x[0] = (log(nu_star) - log(nu_star_MIN[iRAX]))
		/(log(nu_star_MAX[iRAX])-log(nu_star_MIN[iRAX]))*(double)1.8-(double)0.9;// nu*
	x[1]=tanh(pow(GGG,((double)1.0/(double)3.0)))*(double)1.6-(double)0.8;
	x[2]=rho*(double)1.8-(double)0.9;
#ifdef use_NNWwithBETA
	x[3]=beta *(double)8.0 - 0.8;
#endif
	return 0;
}

int CALCULATEbyNNW(struct Weightvalue *pWV, struct NNWvalue *pNNWV, double RAX, double beta)
{
	double net_input;
	int i,j;
	double x_tmp[NUM_INPUT+1]={0}; /* Scaling関数に飛ばすためのダミー変数  */
	NNWVALUE NNWV_tmp;
	
	pNNWV->xunit[NUM_INPUT] = (double)1.0;
	pNNWV->hiddenunit[NUM_HIDDEN] = (double)1.0;

	for ( j = 0; j < NUM_HIDDEN; j++)
	{
		net_input = 0;
		for ( i= 0; i < NUM_INPUT+1; i++)
		{
			net_input = net_input + pWV->w1[i][j] * pNNWV->xunit[i];
		}
		pNNWV->hiddenunit[j] = tanh((double)net_input);
	}
	/* calculation of value of outputs. */
	for ( i = 0; i < NUM_OUTPUT; i++ )
	{
		net_input = 0;
		for ( j = 0; j < NUM_HIDDEN+1; j++)
		{
			net_input = net_input + pWV->w2[j][i] * pNNWV->hiddenunit[j];
		}
		pNNWV->yunit[i] = tanh((double)net_input);
	}

	
	if(NUM_INPUT>1){
		NNWV_tmp.xunit[NUM_INPUT] = (double)1.0;
		NNWV_tmp.hiddenunit[NUM_HIDDEN] = (double)1.0;

		for(i=0;i<NUM_INPUT;i++){
			NNWV_tmp.xunit[i]=pNNWV->xunit[i];
		}
		ScalingInput(x_tmp, 1, 0, 0.5, RAX, beta);
		NNWV_tmp.xunit[1]=x_tmp[1];
		//
		for ( j = 0; j < NUM_HIDDEN; j++)
		{
			net_input = 0;
			for ( i= 0; i < NUM_INPUT+1; i++)
			{
				net_input = net_input + pWV->w1[i][j] * NNWV_tmp.xunit[i];
			}
			NNWV_tmp.hiddenunit[j] = tanh((double)net_input);
		}
		/* calculation of value of outputs.*/
		for ( i = 0; i < NUM_OUTPUT; i++ )
		{
			net_input = 0;
			for ( j = 0; j < NUM_HIDDEN+1; j++)
			{
				net_input = net_input + pWV->w2[j][i] * NNWV_tmp.hiddenunit[j];
			}
			NNWV_tmp.yunit[i] = tanh((double)net_input);
		}
		if(pNNWV->yunit[0] > NNWV_tmp.yunit[0]) pNNWV->yunit[0] = NNWV_tmp.yunit[0];
	}

	return 0;
}
double ERRORFUNCTION(double v_vth)
{
	double PI=	(3.141592653589793e0);
	double herf, herf1;
	herf1 = (DERF(v_vth)-(v_vth*((double)2.0/sqrt(PI))*exp(-(double)1.*(v_vth*v_vth))))
		/((double)2.*(v_vth*v_vth));
	herf=DERF(v_vth)-herf1;
	return(herf); 
}
double DERF(double X)
{
	double PH=1.772453850905516e0;
	double CMC=(double)2.0 / (PH);
	double H=0.5e0;
	double CXC=(double)2.0 * (H)/3.141592653589793e0;
	double CXI=(double)4.0/(CXC);
	
	double V;
	int I;
	double Y;
	int NM = 5;
	int NX = 13;
	int NA = 5;
	double CM[6];
	double CX[14];
	double CQ[14];
	double CA[6];
	double derf;
	double XV;

/*/		DATA CM */
	CM[0] =0.1000000000000000e+01;
	CM[1] =-0.3333333333333333e+00;
	CM[2] =0.1000000000000000e+00;
	CM[3] =-0.2380952380952381e-01;
	CM[4] =0.4629629629629630e-02;
	CM[5] =-0.7575757575757575e-03;
/*/      DATA CX */
	CX[1] = 0.7788007830714048e+00;
    CX[2] =   0.3678794411714423e+00;
	CX[3] =   0.1053992245618643e+00;
    CX[4] =   0.1831563888873418e-01;
    CX[5] =   0.1930454136227709e-02;
    CX[6] =   0.1234098040866796e-03;
    CX[7] =   0.4785117392129009e-05;
    CX[8] =   0.1125351747192591e-06;
    CX[9] =   0.1605228055185612e-08;
    CX[10] =   0.1388794386496402e-10;
    CX[11] =   0.7287724095819692e-13;
    CX[12] =   0.2319522830243569e-15;
    CX[13] =   0.4477732441718302e-18;
/*/  DATA CQ */
	CQ[1] =   0.2500000000000000e+00;
	CQ[2]=   0.1000000000000000e+01;
	CQ[3] =   0.2250000000000000e+01;
	CQ[4] =   0.4000000000000000e+01;
	CQ[5] =   0.6250000000000000e+01;
	CQ[6] =   0.9000000000000000e+01;
	CQ[7] =   0.1225000000000000e+02;
	CQ[8] =   0.1600000000000000e+02;
	CQ[9] =   0.2025000000000000e+02;
	CQ[10] =   0.2500000000000000e+02;
	CQ[11] =   0.3025000000000000e+02;
	CQ[12] =   0.3600000000000000e+02;
	CQ[13] =   0.4225000000000000e+02;
/*/      DATA CA */	
	CA[0] =   0.1000000000000000e+01;
	CA[1] =  -0.1000000000000000e+01;
	CA[2] =   0.3000000000000000e+01;
	CA[3] =  -0.1500000000000000e+02;
	CA[4] =   0.1050000000000000e+03;
	CA[5] =  -0.9449999999999999e+03;

	if(X >= 0)
	{
		XV=X;
	}
    else
	{
		XV=(double)-1.0*X;
	}

	if(XV <= 0.1e0)
	{
		Y = XV*XV;
        V = CM[NM];
        for(I=NM-1;I>-1;I=I-1)
		{
			V = CM[I] + Y * V;
		}
        derf = CMC * XV * V;
		return (derf);
	}
	else if(XV <= 100.0e0)
	{
        Y = XV*XV;
        V = (double)1.0 / ((double)2.0*Y);
        for(I = 1; I<NX+1;I++)
		{
			V = V + CX[I] / (CQ[I] + Y);
		}
        V = CXC * XV * exp((double)-1.0*Y) * V;
        if(XV < 6.0e0) 
		{
			V = V - ((double)2.0/(exp(CXI * XV)-(double)1.0));
		}
        derf = (double)1.0 - V;
	}
	else
	{
		Y = 2 * XV*XV;
		V = CA[NA];
		for (I = NA- 1;I>-1;I=I-1)
		{
			V = CA[I] + Y * V;
		}
		V = exp(-(XV*XV)) / (PH * XV) * V;
		derf = (double)1.0 - V;
	}

    if(X >= 0)
	{
		return(derf);
	}
	else
	{
		derf = -derf;
		return(derf);
	}
	return (derf);
}
void NNWcheck_DCOM_ex(double RAX, double beta)
{
#define NU_G_FILE			"NU_G_R375.dat"
#define NU_RHO_FILE			"NU_RHO_R375.dat"
#define G_RHO_FILE			"G_RHO_R375.dat"
#ifdef use_NNWwithBETA
#define NU_BETA_FILE		"NU_BETA_R375.dat"
#endif

	int i,j,k;
	FILE *NU_G;
	FILE *NU_RHO;
	FILE *G_RHO;
	int iRAX;
	int cnt_num=50;
	double G_const=1e-2;
	double NU_const=1e0;
	double RHO_const=0.5;
	double DATAplot[100][100][5];
	double dNU, dG, dRHO;
	double GGGMIN = 1e-8;
	double RHOMAX = 1.0;
	double RHOMIN = 0.0;
#ifdef use_NNWwithBETA
	FILE *NU_BETA;
	double beta_const=0.00;
	double dbeta;
	double betaMAX=0.05;
	double betaMIN=0.00;
#endif

	printf("NNW check start\nRAX %e: beta %e\n", RAX, beta);
	NU_G=fopen(NU_G_FILE, "w");
	NU_RHO=fopen(NU_RHO_FILE, "w");
	G_RHO=fopen(G_RHO_FILE, "w");
#ifdef use_NNWwithBETA
	NU_BETA=fopen(NU_BETA_FILE, "w");
#endif
	iRAX=(int)(RAX * 100.0);
	switch(iRAX){
	case 345:if(beta==0.0){iRAX = 0;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 350:if(beta==0.0){iRAX = 1;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 353:if(beta==0.0){iRAX = 2;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 360:if(beta==0.0){iRAX = 3;break;}else 
			 if(beta==0.01){iRAX = 6;break;}else 
			 if(beta==0.005){iRAX = 9;break;}else 
			 if(beta==0.02){iRAX = 7;break;}else 
			 if(beta==0.03){iRAX = 8;break;}else 
			 {printf("(3)An improper value was substituted for beta. \n");exit(1);}
	case 375:/*if(beta==0.0){iRAX = 4;break;}else*/
			 if(beta>=0.00 && beta<=0.05){iRAX = 15;break;}else
			 {printf("An improper value was substituted for beta. \n");exit(1);}
	case 390:if(beta==0.0){iRAX = 5;break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	default :printf("An improper value was substituted for RAX. \n");exit(1);
	}

	printf("GGGMAX %e\n", GGGMAX[iRAX]);
	dNU=(log10(nu_star_MAX[iRAX])-log10(nu_star_MIN[iRAX]))/cnt_num;
	dG=(log10(GGGMAX[iRAX])-log10(GGGMIN))/cnt_num;
	dRHO=(RHOMAX-RHOMIN)/cnt_num;
#ifdef use_NNWwithBETA
	dbeta=(betaMAX-betaMIN)/cnt_num;
#endif

	/* NU_RHO */
	for(i=0;i<cnt_num;i++){
		for(j=1;j<cnt_num;j++){
			DATAplot[i][j][0]=pow(10.0, log10(nu_star_MIN[iRAX])+i*dNU);
			DATAplot[i][j][1]=RHOMIN+dRHO*j;
			DATAplot[i][j][2]=G_const;
#ifdef use_NNWwithBETA
			DATAplot[i][j][3]=beta;
#endif
			DATAplot[i][j][NUM_INPUT]=CALCULATE_DIFFUSION_COEFFICIENT(DATAplot[i][j][0], DATAplot[i][j][1], DATAplot[i][j][2], RAX, beta);
		}
	}

	for(i=0;i<cnt_num;i++){
		for(j=1;j<cnt_num;j++){
			for(k=0;k<NUM_INPUT+1;k++)fprintf(NU_RHO, "%e\t",DATAplot[i][j][k]);
			fprintf(NU_RHO, "\n");
		}
		fprintf(NU_RHO, "\n");
	}

	/* NU_G */
	for(i=0;i<cnt_num;i++){
		for(j=0;j<cnt_num;j++){
			DATAplot[i][j][0]=pow(10.0, (log10(nu_star_MIN[iRAX])+i*dNU));
			DATAplot[i][j][1]=RHO_const;
			DATAplot[i][j][2]=pow(10.0, log10(GGGMIN)+j*dG);
#ifdef use_NNWwithBETA
			DATAplot[i][j][3]=beta;
#endif
			DATAplot[i][j][NUM_INPUT]=CALCULATE_DIFFUSION_COEFFICIENT(DATAplot[i][j][0], DATAplot[i][j][1], DATAplot[i][j][2], RAX, beta);
		}
	}
	for(i=0;i<cnt_num;i++){
		for(j=0;j<cnt_num;j++){
			fprintf(NU_G, "%e\t%e\t%e\t%e\n",DATAplot[i][j][0],DATAplot[i][j][1],DATAplot[i][j][2],DATAplot[i][j][NUM_INPUT]);
		}
		fprintf(NU_G, "\n");
	}

	/* G_RHO */
	for(i=1;i<cnt_num;i++){
		for(j=0;j<cnt_num;j++){
			DATAplot[i][j][0]=NU_const;
			DATAplot[i][j][1]=RHOMIN+dRHO*i;
			DATAplot[i][j][2]=pow(10.0, log10(GGGMIN)+j*dG);
#ifdef use_NNWwithBETA
			DATAplot[i][j][3]=beta;
#endif
			DATAplot[i][j][NUM_INPUT]=CALCULATE_DIFFUSION_COEFFICIENT(DATAplot[i][j][0], DATAplot[i][j][1], DATAplot[i][j][2], RAX, beta);
		}
	}
	for(i=1;i<cnt_num;i++){
		for(j=0;j<cnt_num;j++){
			fprintf(G_RHO, "%e\t%e\t%e\t%e\n",DATAplot[i][j][0],DATAplot[i][j][1],DATAplot[i][j][2],DATAplot[i][j][NUM_INPUT]);
		}
		fprintf(G_RHO, "\n");
	}
	fclose(NU_G);
	fclose(G_RHO);
	fclose(NU_RHO);
	printf("output end of NNW check\n");
	exit(0);
}
int NNWSETELMENT(){
/* Setting of NNW data of each magnetic field
// iRAX NUMBER (beta = 0%)
//      RAX345-> 0
//      RAX350-> 1
//		RAX353-> 2
//		RAX360-> 3
//		RAX375-> 4
//		RAX390-> 5
// iRAX NUMBER (finite beta)  
//      RAX360 beta = 0.5% -> 9
//      RAX360 beta = 1% -> 6
//      RAX360 beta = 2% -> 7
//      RAX360 beta = 3% -> 8
*/
/////////////////   RAX 345 data ////////////////
WV[0].omega[0] = 	-1.60445E+00	;
WV[0].omega[1] = 	-2.27848E+00	;
WV[0].omega[2] = 	-3.00908E+00	;
WV[0].omega[3] = 	-2.06735E+00	;
WV[0].omega[4] = 	2.10436E+00	;
WV[0].omega[5] = 	-2.56054E+00	;
WV[0].omega[6] = 	-3.48413E+00	;
WV[0].omega[7] = 	-3.42630E+00	;
WV[0].omega[8] = 	-4.04821E+00	;
WV[0].omega[9] = 	1.56369E+00	;
WV[0].omega[10] =	-4.80120E+00	;
WV[0].omega[11] =	-2.13940E+00	;
WV[0].omega[12] =	-3.23820E+00	;
WV[0].omega[13] =	-5.34890E+00	;
WV[0].omega[14] =	7.95560E-01	;
WV[0].omega[15] =	-7.56682E+00	;
WV[0].omega[16] =	1.62098E+00	;
WV[0].omega[17] =	-1.91408E+01	;
WV[0].omega[18] =	-4.76747E-02	;
WV[0].omega[19] =	-3.23878E-01	;
WV[0].omega[20] =	-9.98009E+00	;
WV[0].omega[21] =	1.93451E+00	;
WV[0].omega[22] =	1.82593E+00	;
WV[0].omega[23] =	2.19065E-01	;
WV[0].omega[24] =	5.20978E+00	;
WV[0].omega[25] =	8.05802E-01	;
WV[0].omega[26] =	1.47081E+00	;
WV[0].omega[27] =	6.23130E-03	;
WV[0].omega[28] =	3.24522E+00	;
WV[0].omega[29] =	-5.62248E+00	;
WV[0].omega[30] =	2.20214E-01	;
WV[0].omega[31] =	1.15158E+00	;
WV[0].omega[32] =	8.36018E-01	;
WV[0].omega[33] =	7.18371E-02	;
WV[0].omega[34] =	1.34044E+00	;
WV[0].omega[35] =	3.46518E-01	;
WV[0].omega[36] =	-2.45781E+00	;
WV[0].omega[37] =	-2.33139E+00	;
WV[0].omega[38] =	-9.47805E-02	;
WV[0].omega[39] =	-1.47189E-01	;
WV[0].omega[40] =	-3.50581E-01	;
WV[0].omega[41] =	1.36029E+00	;
WV[0].omega[42] =	1.47177E+00	;
WV[0].omega[43] =	3.97300E-01	;
WV[0].omega[44] =	6.59612E+00	;
WV[0].omega[45] =	-7.52496E+00	;
WV[0].omega[46] =	-1.13507E+00	;
WV[0].omega[47] =	-1.15590E+01	;
WV[0].omega[48] =	4.76596E-02	;
WV[0].omega[49] =	8.60942E-01	;
WV[0].omega[50] =	-9.97976E+00	;
WV[0].omega[51] =	-1.74335E+00	;
WV[0].omega[52] =	-1.71011E+00	;
WV[0].omega[53] =	8.99666E-02	;
WV[0].omega[54] =	4.39308E+00	;
WV[0].omega[55] =	-1.01634E+00	;
WV[0].omega[56] =	-1.18316E+00	;
WV[0].omega[57] =	-7.54882E-01	;
WV[0].omega[58] =	-5.63056E-01	;
WV[0].omega[59] =	-8.21386E+00	;
WV[0].omega[60] =	-4.32640E+00	;
WV[0].omega[61] =	-2.15096E+00	;
WV[0].omega[62] =	-2.54328E-02	;
WV[0].omega[63] =	-1.01060E+00	;
WV[0].omega[64] =	1.51946E-01	;
WV[0].omega[65] =	3.31666E+00	;
WV[0].omega[66] =	-1.99017E+00	;
WV[0].omega[67] =	2.12066E+00	;
WV[0].omega[68] =	4.85115E-01	;
WV[0].omega[69] =	-5.57337E-01	;
WV[0].omega[70] =	3.44796E-01	;
WV[0].omega[71] =	1.80447E+00	;
WV[0].omega[72] =	1.74235E-01	;
WV[0].omega[73] =	-1.99033E-01	;
WV[0].omega[74] =	2.27678E-02	;
WV[0].omega[75] =	-6.20471E-01	;
WV[0].w1[0][0] = 	-1.60445E+00	;
WV[0].w1[0][1] = 	-2.27848E+00	;
WV[0].w1[0][2] = 	-3.00908E+00	;
WV[0].w1[0][3] = 	-2.06735E+00	;
WV[0].w1[0][4] = 	2.10436E+00	;
WV[0].w1[0][5] = 	-2.56054E+00	;
WV[0].w1[0][6] = 	-3.48413E+00	;
WV[0].w1[0][7] = 	-3.42630E+00	;
WV[0].w1[0][8] = 	-4.04821E+00	;
WV[0].w1[0][9] = 	1.56369E+00	;
WV[0].w1[0][10] =	-4.80120E+00	;
WV[0].w1[0][11] =	-2.13940E+00	;
WV[0].w1[0][12] =	-3.23820E+00	;
WV[0].w1[0][13] =	-5.34890E+00	;
WV[0].w1[0][14] =	7.95560E-01	;
WV[0].w1[1][0] = 	-7.56682E+00	;
WV[0].w1[1][1] = 	1.62098E+00	;
WV[0].w1[1][2] = 	-1.91408E+01	;
WV[0].w1[1][3] = 	-4.76747E-02	;
WV[0].w1[1][4] = 	-3.23878E-01	;
WV[0].w1[1][5] = 	-9.98009E+00	;
WV[0].w1[1][6] = 	1.93451E+00	;
WV[0].w1[1][7] = 	1.82593E+00	;
WV[0].w1[1][8] = 	2.19065E-01	;
WV[0].w1[1][9] = 	5.20978E+00	;
WV[0].w1[1][10] =	8.05802E-01	;
WV[0].w1[1][11] =	1.47081E+00	;
WV[0].w1[1][12] =	6.23130E-03	;
WV[0].w1[1][13] =	3.24522E+00	;
WV[0].w1[1][14] =	-5.62248E+00	;
WV[0].w1[2][0] = 	2.20214E-01	;
WV[0].w1[2][1] = 	1.15158E+00	;
WV[0].w1[2][2] = 	8.36018E-01	;
WV[0].w1[2][3] = 	7.18371E-02	;
WV[0].w1[2][4] = 	1.34044E+00	;
WV[0].w1[2][5] = 	3.46518E-01	;
WV[0].w1[2][6] = 	-2.45781E+00	;
WV[0].w1[2][7] = 	-2.33139E+00	;
WV[0].w1[2][8] = 	-9.47805E-02	;
WV[0].w1[2][9] = 	-1.47189E-01	;
WV[0].w1[2][10] =	-3.50581E-01	;
WV[0].w1[2][11] =	1.36029E+00	;
WV[0].w1[2][12] =	1.47177E+00	;
WV[0].w1[2][13] =	3.97300E-01	;
WV[0].w1[2][14] =	6.59612E+00	;
WV[0].w1[3][0] = 	-7.52496E+00	;
WV[0].w1[3][1] = 	-1.13507E+00	;
WV[0].w1[3][2] = 	-1.15590E+01	;
WV[0].w1[3][3] = 	4.76596E-02	;
WV[0].w1[3][4] = 	8.60942E-01	;
WV[0].w1[3][5] = 	-9.97976E+00	;
WV[0].w1[3][6] = 	-1.74335E+00	;
WV[0].w1[3][7] = 	-1.71011E+00	;
WV[0].w1[3][8] = 	8.99666E-02	;
WV[0].w1[3][9] = 	4.39308E+00	;
WV[0].w1[3][10] =	-1.01634E+00	;
WV[0].w1[3][11] =	-1.18316E+00	;
WV[0].w1[3][12] =	-7.54882E-01	;
WV[0].w1[3][13] =	-5.63056E-01	;
WV[0].w1[3][14] =	-8.21386E+00	;
WV[0].w2[0][0] = 	-4.32640E+00	;
WV[0].w2[1][0] = 	-2.15096E+00	;
WV[0].w2[2][0] = 	-2.54328E-02	;
WV[0].w2[3][0] = 	-1.01060E+00	;
WV[0].w2[4][0] = 	1.51946E-01	;
WV[0].w2[5][0] = 	3.31666E+00	;
WV[0].w2[6][0] = 	-1.99017E+00	;
WV[0].w2[7][0] = 	2.12066E+00	;
WV[0].w2[8][0] = 	4.85115E-01	;
WV[0].w2[9][0] = 	-5.57337E-01	;
WV[0].w2[10][0] =	3.44796E-01	;
WV[0].w2[11][0] =	1.80447E+00	;
WV[0].w2[12][0] =	1.74235E-01	;
WV[0].w2[13][0] =	-1.99033E-01	;
WV[0].w2[14][0] =	2.27678E-02	;
WV[0].w2[15][0] =	-6.20471E-01	;
nu_star_MAX[0] = 1.000000e+003;	
nu_star_MIN[0] = 1.000000e-006;	
D_star_MAX[0]  = 1.000000e+004;	
D_star_MIN[0]  = 1.000000e-004;	
GGGMAX[0]      = 1.000000e-001;	

/////////////////   RAX 350 data ////////////////
WV[1].omega[0] = 	-3.08786E+00	;
WV[1].omega[1] = 	-4.71454E+00	;
WV[1].omega[2] = 	-2.73297E+00	;
WV[1].omega[3] = 	2.10422E+00	;
WV[1].omega[4] = 	-2.19632E+00	;
WV[1].omega[5] = 	-1.65395E+00	;
WV[1].omega[6] = 	2.42630E+00	;
WV[1].omega[7] = 	-1.21853E+00	;
WV[1].omega[8] = 	2.33939E+00	;
WV[1].omega[9] = 	9.54569E+00	;
WV[1].omega[10] =	-2.08731E+00	;
WV[1].omega[11] =	8.65086E-01	;
WV[1].omega[12] =	-6.10750E+00	;
WV[1].omega[13] =	3.27668E+00	;
WV[1].omega[14] =	2.90368E+00	;
WV[1].omega[15] =	6.06835E-01	;
WV[1].omega[16] =	3.71016E+00	;
WV[1].omega[17] =	-1.40070E+00	;
WV[1].omega[18] =	-4.00456E-01	;
WV[1].omega[19] =	1.41934E+00	;
WV[1].omega[20] =	-2.84708E-01	;
WV[1].omega[21] =	7.05384E+00	;
WV[1].omega[22] =	1.28175E+00	;
WV[1].omega[23] =	-6.89963E-01	;
WV[1].omega[24] =	6.69946E-02	;
WV[1].omega[25] =	1.52303E-02	;
WV[1].omega[26] =	-6.51106E-01	;
WV[1].omega[27] =	9.59596E+00	;
WV[1].omega[28] =	4.46527E+00	;
WV[1].omega[29] =	2.01750E+00	;
WV[1].omega[30] =	4.44719E-02	;
WV[1].omega[31] =	8.08182E-01	;
WV[1].omega[32] =	5.30346E-01	;
WV[1].omega[33] =	-3.59380E-01	;
WV[1].omega[34] =	7.86589E+00	;
WV[1].omega[35] =	-8.73060E-01	;
WV[1].omega[36] =	-1.39705E-01	;
WV[1].omega[37] =	-1.60627E+00	;
WV[1].omega[38] =	-5.94957E-01	;
WV[1].omega[39] =	6.68244E-01	;
WV[1].omega[40] =	-4.55088E-02	;
WV[1].omega[41] =	8.18950E-01	;
WV[1].omega[42] =	1.03065E+00	;
WV[1].omega[43] =	-2.72381E-01	;
WV[1].omega[44] =	-1.87039E+00	;
WV[1].omega[45] =	-6.56178E-01	;
WV[1].omega[46] =	-8.42167E-01	;
WV[1].omega[47] =	-2.15679E+00	;
WV[1].omega[48] =	7.46450E-01	;
WV[1].omega[49] =	-3.48418E+00	;
WV[1].omega[50] =	-6.45031E-01	;
WV[1].omega[51] =	8.69207E+00	;
WV[1].omega[52] =	-1.28088E+00	;
WV[1].omega[53] =	6.90501E-01	;
WV[1].omega[54] =	-1.61023E-01	;
WV[1].omega[55] =	2.03354E+00	;
WV[1].omega[56] =	2.23605E+00	;
WV[1].omega[57] =	2.22430E+00	;
WV[1].omega[58] =	6.61244E+00	;
WV[1].omega[59] =	3.05211E+00	;
WV[1].omega[60] =	7.62590E-01	;
WV[1].omega[61] =	-1.85913E-01	;
WV[1].omega[62] =	4.93663E-01	;
WV[1].omega[63] =	2.97738E+00	;
WV[1].omega[64] =	1.91762E-02	;
WV[1].omega[65] =	-1.41679E-01	;
WV[1].omega[66] =	6.15744E+00	;
WV[1].omega[67] =	-2.91147E-01	;
WV[1].omega[68] =	-1.60622E+00	;
WV[1].omega[69] =	-3.48102E-02	;
WV[1].omega[70] =	-8.45378E-01	;
WV[1].omega[71] =	-4.34289E+00	;
WV[1].omega[72] =	-1.74437E-02	;
WV[1].omega[73] =	-3.09505E+00	;
WV[1].omega[74] =	6.18897E-02	;
WV[1].omega[75] =	1.30198E+00	;
WV[1].w1[0][0] = 	-3.08786E+00	;
WV[1].w1[0][1] = 	-4.71454E+00	;
WV[1].w1[0][2] = 	-2.73297E+00	;
WV[1].w1[0][3] = 	2.10422E+00	;
WV[1].w1[0][4] = 	-2.19632E+00	;
WV[1].w1[0][5] = 	-1.65395E+00	;
WV[1].w1[0][6] = 	2.42630E+00	;
WV[1].w1[0][7] = 	-1.21853E+00	;
WV[1].w1[0][8] = 	2.33939E+00	;
WV[1].w1[0][9] = 	9.54569E+00	;
WV[1].w1[0][10] =	-2.08731E+00	;
WV[1].w1[0][11] =	8.65086E-01	;
WV[1].w1[0][12] =	-6.10750E+00	;
WV[1].w1[0][13] =	3.27668E+00	;
WV[1].w1[0][14] =	2.90368E+00	;
WV[1].w1[1][0] = 	6.06835E-01	;
WV[1].w1[1][1] = 	3.71016E+00	;
WV[1].w1[1][2] = 	-1.40070E+00	;
WV[1].w1[1][3] = 	-4.00456E-01	;
WV[1].w1[1][4] = 	1.41934E+00	;
WV[1].w1[1][5] = 	-2.84708E-01	;
WV[1].w1[1][6] = 	7.05384E+00	;
WV[1].w1[1][7] = 	1.28175E+00	;
WV[1].w1[1][8] = 	-6.89963E-01	;
WV[1].w1[1][9] = 	6.69946E-02	;
WV[1].w1[1][10] =	1.52303E-02	;
WV[1].w1[1][11] =	-6.51106E-01	;
WV[1].w1[1][12] =	9.59596E+00	;
WV[1].w1[1][13] =	4.46527E+00	;
WV[1].w1[1][14] =	2.01750E+00	;
WV[1].w1[2][0] = 	4.44719E-02	;
WV[1].w1[2][1] = 	8.08182E-01	;
WV[1].w1[2][2] = 	5.30346E-01	;
WV[1].w1[2][3] = 	-3.59380E-01	;
WV[1].w1[2][4] = 	7.86589E+00	;
WV[1].w1[2][5] = 	-8.73060E-01	;
WV[1].w1[2][6] = 	-1.39705E-01	;
WV[1].w1[2][7] = 	-1.60627E+00	;
WV[1].w1[2][8] = 	-5.94957E-01	;
WV[1].w1[2][9] = 	6.68244E-01	;
WV[1].w1[2][10] =	-4.55088E-02	;
WV[1].w1[2][11] =	8.18950E-01	;
WV[1].w1[2][12] =	1.03065E+00	;
WV[1].w1[2][13] =	-2.72381E-01	;
WV[1].w1[2][14] =	-1.87039E+00	;
WV[1].w1[3][0] = 	-6.56178E-01	;
WV[1].w1[3][1] = 	-8.42167E-01	;
WV[1].w1[3][2] = 	-2.15679E+00	;
WV[1].w1[3][3] = 	7.46450E-01	;
WV[1].w1[3][4] = 	-3.48418E+00	;
WV[1].w1[3][5] = 	-6.45031E-01	;
WV[1].w1[3][6] = 	8.69207E+00	;
WV[1].w1[3][7] = 	-1.28088E+00	;
WV[1].w1[3][8] = 	6.90501E-01	;
WV[1].w1[3][9] = 	-1.61023E-01	;
WV[1].w1[3][10] =	2.03354E+00	;
WV[1].w1[3][11] =	2.23605E+00	;
WV[1].w1[3][12] =	2.22430E+00	;
WV[1].w1[3][13] =	6.61244E+00	;
WV[1].w1[3][14] =	3.05211E+00	;
WV[1].w2[0][0] = 	7.62590E-01	;
WV[1].w2[1][0] = 	-1.85913E-01	;
WV[1].w2[2][0] = 	4.93663E-01	;
WV[1].w2[3][0] = 	2.97738E+00	;
WV[1].w2[4][0] = 	1.91762E-02	;
WV[1].w2[5][0] = 	-1.41679E-01	;
WV[1].w2[6][0] = 	6.15744E+00	;
WV[1].w2[7][0] = 	-2.91147E-01	;
WV[1].w2[8][0] = 	-1.60622E+00	;
WV[1].w2[9][0] = 	-3.48102E-02	;
WV[1].w2[10][0] =	-8.45378E-01	;
WV[1].w2[11][0] =	-4.34289E+00	;
WV[1].w2[12][0] =	-1.74437E-02	;
WV[1].w2[13][0] =	-3.09505E+00	;
WV[1].w2[14][0] =	6.18897E-02	;
WV[1].w2[15][0] =	1.30198E+00	;
nu_star_MAX[1] = 1.000000e+003;	
nu_star_MIN[1] = 1.000000e-006;	
D_star_MAX[1]  = 1.000000e+005;	
D_star_MIN[1]  = 1.000000e-004;	
GGGMAX[1]      = 1.000000e-001;	

/////////////////   RAX 353 data ////////////////
WV[2].omega[0] = -1.356671e+000;
WV[2].omega[1] = -1.561201e-001;
WV[2].omega[2] = -1.807516e+000;
WV[2].omega[3] = 8.128963e+000;
WV[2].omega[4] = 4.186174e+000;
WV[2].omega[5] = 2.830330e+000;
WV[2].omega[6] = 2.132897e+000;
WV[2].omega[7] = 4.113541e-001;
WV[2].omega[8] = -6.891069e+000;
WV[2].omega[9] = -4.711259e+000;
WV[2].omega[10] =-3.288491e-001;
WV[2].omega[11] =5.195396e+000;
WV[2].omega[12] =-1.808857e+000;
WV[2].omega[13] =1.421846e+000;
WV[2].omega[14] =-1.915895e+000;
WV[2].omega[15] =2.499396e-002;
WV[2].omega[16] =-1.957575e+000;
WV[2].omega[17] =2.265711e+000;
WV[2].omega[18] =-3.151872e-001;
WV[2].omega[19] =-1.353183e-001;
WV[2].omega[20] =4.093474e-002;
WV[2].omega[21] =-2.205593e+000;
WV[2].omega[22] =-2.398775e-001;
WV[2].omega[23] =-5.838076e-002;
WV[2].omega[24] =6.827993e-003;
WV[2].omega[25] =-1.590913e+000;
WV[2].omega[26] =5.298194e-002;
WV[2].omega[27] =1.495882e-002;
WV[2].omega[28] =-1.303338e-002;
WV[2].omega[29] =6.257244e-002;
WV[2].omega[30] =2.423047e+000;
WV[2].omega[31] =-5.470278e-002;
WV[2].omega[32] =2.904401e-001;
WV[2].omega[33] =1.158586e-001;
WV[2].omega[34] =-2.842489e-001;
WV[2].omega[35] =5.867862e-001;
WV[2].omega[36] =-3.217529e-001;
WV[2].omega[37] =-3.516647e-001;
WV[2].omega[38] =1.952571e+000;
WV[2].omega[39] =-1.673365e-001;
WV[2].omega[40] =-2.669276e-002;
WV[2].omega[41] =3.620579e-001;
WV[2].omega[42] =-1.178359e+000;
WV[2].omega[43] =-2.759827e+000;
WV[2].omega[44] =-1.645364e+000;
WV[2].omega[45] =3.710146e-001;
WV[2].omega[46] =-1.279781e+000;
WV[2].omega[47] =1.806111e-001;
WV[2].omega[48] =3.666050e+000;
WV[2].omega[49] =-6.089531e-001;
WV[2].omega[50] =-1.117915e+000;
WV[2].omega[51] =-1.853228e-001;
WV[2].omega[52] =-1.493609e-001;
WV[2].omega[53] =-2.174618e+000;
WV[2].omega[54] =6.438117e-001;
WV[2].omega[55] =-9.645149e-001;
WV[2].omega[56] =6.466592e-001;
WV[2].omega[57] =4.059488e-001;
WV[2].omega[58] =-3.985659e-001;
WV[2].omega[59] =3.828265e-001;
WV[2].omega[60] =-1.006650e+000;
WV[2].omega[61] =1.816740e+000;
WV[2].omega[62] =-2.880433e+000;
WV[2].omega[63] =6.243385e-002;
WV[2].omega[64] =-4.696891e-001;
WV[2].omega[65] =6.717203e-001;
WV[2].omega[66] =-2.499011e+000;
WV[2].omega[67] =-1.013723e+000;
WV[2].omega[68] =6.893917e-002;
WV[2].omega[69] =-3.578033e-001;
WV[2].omega[70] =-2.134360e+000;
WV[2].omega[71] =1.914079e-001;
WV[2].omega[72] =1.177335e+000;
WV[2].omega[73] =-7.881831e-001;
WV[2].omega[74] =-5.662642e-001;
WV[2].omega[75] =-1.490022e-001;
WV[2].w1[0][0] = -1.356671e+000;
WV[2].w1[0][1] = -1.561201e-001;
WV[2].w1[0][2] = -1.807516e+000;
WV[2].w1[0][3] = 8.128963e+000;
WV[2].w1[0][4] = 4.186174e+000;
WV[2].w1[0][5] = 2.830330e+000;
WV[2].w1[0][6] = 2.132897e+000;
WV[2].w1[0][7] = 4.113541e-001;
WV[2].w1[0][8] = -6.891069e+000;
WV[2].w1[0][9] = -4.711259e+000;
WV[2].w1[0][10] =-3.288491e-001;
WV[2].w1[0][11] =5.195396e+000;
WV[2].w1[0][12] =-1.808857e+000;
WV[2].w1[0][13] =1.421846e+000;
WV[2].w1[0][14] =-1.915895e+000;
WV[2].w1[1][0] = 2.499396e-002;
WV[2].w1[1][1] = -1.957575e+000;
WV[2].w1[1][2] = 2.265711e+000;
WV[2].w1[1][3] = -3.151872e-001;
WV[2].w1[1][4] = -1.353183e-001;
WV[2].w1[1][5] = 4.093474e-002;
WV[2].w1[1][6] = -2.205593e+000;
WV[2].w1[1][7] = -2.398775e-001;
WV[2].w1[1][8] = -5.838076e-002;
WV[2].w1[1][9] = 6.827993e-003;
WV[2].w1[1][10] =-1.590913e+000;
WV[2].w1[1][11] =5.298194e-002;
WV[2].w1[1][12] =1.495882e-002;
WV[2].w1[1][13] =-1.303338e-002;
WV[2].w1[1][14] =6.257244e-002;
WV[2].w1[2][0] = 2.423047e+000;
WV[2].w1[2][1] = -5.470278e-002;
WV[2].w1[2][2] = 2.904401e-001;
WV[2].w1[2][3] = 1.158586e-001;
WV[2].w1[2][4] = -2.842489e-001;
WV[2].w1[2][5] = 5.867862e-001;
WV[2].w1[2][6] = -3.217529e-001;
WV[2].w1[2][7] = -3.516647e-001;
WV[2].w1[2][8] = 1.952571e+000;
WV[2].w1[2][9] = -1.673365e-001;
WV[2].w1[2][10] =-2.669276e-002;
WV[2].w1[2][11] =3.620579e-001;
WV[2].w1[2][12] =-1.178359e+000;
WV[2].w1[2][13] =-2.759827e+000;
WV[2].w1[2][14] =-1.645364e+000;
WV[2].w1[3][0] = 3.710146e-001;
WV[2].w1[3][1] = -1.279781e+000;
WV[2].w1[3][2] = 1.806111e-001;
WV[2].w1[3][3] = 3.666050e+000;
WV[2].w1[3][4] = -6.089531e-001;
WV[2].w1[3][5] = -1.117915e+000;
WV[2].w1[3][6] = -1.853228e-001;
WV[2].w1[3][7] = -1.493609e-001;
WV[2].w1[3][8] = -2.174618e+000;
WV[2].w1[3][9] = 6.438117e-001;
WV[2].w1[3][10] =-9.645149e-001;
WV[2].w1[3][11] =6.466592e-001;
WV[2].w1[3][12] =4.059488e-001;
WV[2].w1[3][13] =-3.985659e-001;
WV[2].w1[3][14] =3.828265e-001;
WV[2].w2[0][0] = -1.006650e+000;
WV[2].w2[1][0] = 1.816740e+000;
WV[2].w2[2][0] = -2.880433e+000;
WV[2].w2[3][0] = 6.243385e-002;
WV[2].w2[4][0] = -4.696891e-001;
WV[2].w2[5][0] = 6.717203e-001;
WV[2].w2[6][0] = -2.499011e+000;
WV[2].w2[7][0] = -1.013723e+000;
WV[2].w2[8][0] = 6.893917e-002;
WV[2].w2[9][0] = -3.578033e-001;
WV[2].w2[10][0] =-2.134360e+000;
WV[2].w2[11][0] =1.914079e-001;
WV[2].w2[12][0] =1.177335e+000;
WV[2].w2[13][0] =-7.881831e-001;
WV[2].w2[14][0] =-5.662642e-001;
WV[2].w2[15][0] =-1.490022e-001;
nu_star_MAX[2] = 1.000000e+003;
nu_star_MIN[2] = 1.000000e-006;
D_star_MAX[2]  = 1.000000e+005;
D_star_MIN[2]  = 1.000000e-004;
GGGMAX[2]      = 0.800000e-001;

/////////////////   RAX 360 data ////////////////

#ifdef NEW_NNW_20121019
    WV[3].omega[0]	=	-2.09509E+00	;
    WV[3].omega[1]	=	-5.55139E+00	;
    WV[3].omega[2]	=	-8.03140E+00	;
    WV[3].omega[3]	=	-1.33101E+01	;
    WV[3].omega[4]	=	2.08584E+00	;
    WV[3].omega[5]	=	-2.40708E+00	;
    WV[3].omega[6]	=	-5.20571E+00	;
    WV[3].omega[7]	=	-1.55396E+00	;
    WV[3].omega[8]	=	8.95843E+00	;
    WV[3].omega[9]	=	-9.58097E+00	;
    WV[3].omega[10]	=	-1.70512E+00	;
    WV[3].omega[11]	=	-9.18579E+00	;
    WV[3].omega[12]	=	-9.27721E+00	;
    WV[3].omega[13]	=	5.63908E+00	;
    WV[3].omega[14]	=	1.92097E+00	;
    WV[3].omega[15]	=	1.05606E+00	;
    WV[3].omega[16]	=	9.08925E-01	;
    WV[3].omega[17]	=	1.98181E-01	;
    WV[3].omega[18]	=	-1.80701E-01	;
    WV[3].omega[19]	=	-8.96149E-01	;
    WV[3].omega[20]	=	-2.95224E-01	;
    WV[3].omega[21]	=	6.70643E-01	;
    WV[3].omega[22]	=	-3.25353E-01	;
    WV[3].omega[23]	=	1.23043E-01	;
    WV[3].omega[24]	=	2.75193E-01	;
    WV[3].omega[25]	=	-4.51169E-01	;
    WV[3].omega[26]	=	5.34783E-02	;
    WV[3].omega[27]	=	-2.50249E+01	;
    WV[3].omega[28]	=	1.19552E+01	;
    WV[3].omega[29]	=	-3.85880E-01	;
    WV[3].omega[30]	=	7.68857E-01	;
    WV[3].omega[31]	=	3.85240E-01	;
    WV[3].omega[32]	=	-5.34414E-01	;
    WV[3].omega[33]	=	-5.27117E-02	;
    WV[3].omega[34]	=	-7.40683E-01	;
    WV[3].omega[35]	=	3.85173E-01	;
    WV[3].omega[36]	=	3.70093E-01	;
    WV[3].omega[37]	=	-1.78522E-01	;
    WV[3].omega[38]	=	-3.04144E-01	;
    WV[3].omega[39]	=	-6.74837E-01	;
    WV[3].omega[40]	=	-1.58832E-01	;
    WV[3].omega[41]	=	-1.90543E-01	;
    WV[3].omega[42]	=	8.89136E-01	;
    WV[3].omega[43]	=	-3.69381E-01	;
    WV[3].omega[44]	=	-7.06841E-01	;
    WV[3].omega[45]	=	-2.14645E-01	;
    WV[3].omega[46]	=	-3.13397E-01	;
    WV[3].omega[47]	=	-4.80665E+00	;
    WV[3].omega[48]	=	2.26589E+00	;
    WV[3].omega[49]	=	2.53120E-01	;
    WV[3].omega[50]	=	-4.59605E-01	;
    WV[3].omega[51]	=	-3.22872E-01	;
    WV[3].omega[52]	=	-1.89153E+00	;
    WV[3].omega[53]	=	-1.26296E+00	;
    WV[3].omega[54]	=	-5.06942E+00	;
    WV[3].omega[55]	=	-1.61459E+00	;
    WV[3].omega[56]	=	-4.00539E+00	;
    WV[3].omega[57]	=	-2.10757E+01	;
    WV[3].omega[58]	=	1.15692E+01	;
    WV[3].omega[59]	=	4.06940E-01	;
    WV[3].omega[60]	=	1.04884E+01	;
    WV[3].omega[61]	=	-3.32807E+00	;
    WV[3].omega[62]	=	6.11381E+01	;
    WV[3].omega[63]	=	-6.89613E-02	;
    WV[3].omega[64]	=	1.54260E+01	;
    WV[3].omega[65]	=	-2.34379E+00	;
    WV[3].omega[66]	=	4.51307E+00	;
    WV[3].omega[67]	=	-3.88168E+01	;
    WV[3].omega[68]	=	-2.09305E-01	;
    WV[3].omega[69]	=	-2.28957E+01	;
    WV[3].omega[70]	=	1.83818E+01	;
    WV[3].omega[71]	=	-1.20588E+00	;
    WV[3].omega[72]	=	-1.07955E-01	;
    WV[3].omega[73]	=	-5.73910E-01	;
    WV[3].omega[74]	=	-5.30859E+00	;
    WV[3].omega[75]	=	1.70340E+01	;
    WV[3].w1[0][0]	=	-2.09509E+00	;
    WV[3].w1[0][1]	=	-5.55139E+00	;
    WV[3].w1[0][2]	=	-8.03140E+00	;
    WV[3].w1[0][3]	=	-1.33101E+01	;
    WV[3].w1[0][4]	=	2.08584E+00	;
    WV[3].w1[0][5]	=	-2.40708E+00	;
    WV[3].w1[0][6]	=	-5.20571E+00	;
    WV[3].w1[0][7]	=	-1.55396E+00	;
    WV[3].w1[0][8]	=	8.95843E+00	;
    WV[3].w1[0][9]	=	-9.58097E+00	;
    WV[3].w1[0][10]	=	-1.70512E+00	;
    WV[3].w1[0][11]	=	-9.18579E+00	;
    WV[3].w1[0][12]	=	-9.27721E+00	;
    WV[3].w1[0][13]	=	5.63908E+00	;
    WV[3].w1[0][14]	=	1.92097E+00	;
    WV[3].w1[1][0]	=	1.05606E+00	;
    WV[3].w1[1][1]	=	9.08925E-01	;
    WV[3].w1[1][2]	=	1.98181E-01	;
    WV[3].w1[1][3]	=	-1.80701E-01	;
    WV[3].w1[1][4]	=	-8.96149E-01	;
    WV[3].w1[1][5]	=	-2.95224E-01	;
    WV[3].w1[1][6]	=	6.70643E-01	;
    WV[3].w1[1][7]	=	-3.25353E-01	;
    WV[3].w1[1][8]	=	1.23043E-01	;
    WV[3].w1[1][9]	=	2.75193E-01	;
    WV[3].w1[1][10]	=	-4.51169E-01	;
    WV[3].w1[1][11]	=	5.34783E-02	;
    WV[3].w1[1][12]	=	-2.50249E+01	;
    WV[3].w1[1][13]	=	1.19552E+01	;
    WV[3].w1[1][14]	=	-3.85880E-01	;
    WV[3].w1[2][0]	=	7.68857E-01	;
    WV[3].w1[2][1]	=	3.85240E-01	;
    WV[3].w1[2][2]	=	-5.34414E-01	;
    WV[3].w1[2][3]	=	-5.27117E-02	;
    WV[3].w1[2][4]	=	-7.40683E-01	;
    WV[3].w1[2][5]	=	3.85173E-01	;
    WV[3].w1[2][6]	=	3.70093E-01	;
    WV[3].w1[2][7]	=	-1.78522E-01	;
    WV[3].w1[2][8]	=	-3.04144E-01	;
    WV[3].w1[2][9]	=	-6.74837E-01	;
    WV[3].w1[2][10]	=	-1.58832E-01	;
    WV[3].w1[2][11]	=	-1.90543E-01	;
    WV[3].w1[2][12]	=	8.89136E-01	;
    WV[3].w1[2][13]	=	-3.69381E-01	;
    WV[3].w1[2][14]	=	-7.06841E-01	;
    WV[3].w1[3][0]	=	-2.14645E-01	;
    WV[3].w1[3][1]	=	-3.13397E-01	;
    WV[3].w1[3][2]	=	-4.80665E+00	;
    WV[3].w1[3][3]	=	2.26589E+00	;
    WV[3].w1[3][4]	=	2.53120E-01	;
    WV[3].w1[3][5]	=	-4.59605E-01	;
    WV[3].w1[3][6]	=	-3.22872E-01	;
    WV[3].w1[3][7]	=	-1.89153E+00	;
    WV[3].w1[3][8]	=	-1.26296E+00	;
    WV[3].w1[3][9]	=	-5.06942E+00	;
    WV[3].w1[3][10]	=	-1.61459E+00	;
    WV[3].w1[3][11]	=	-4.00539E+00	;
    WV[3].w1[3][12]	=	-2.10757E+01	;
    WV[3].w1[3][13]	=	1.15692E+01	;
    WV[3].w1[3][14]	=	4.06940E-01	;
    WV[3].w2[0][0]	=	1.04884E+01	;
    WV[3].w2[1][0]	=	-3.32807E+00	;
    WV[3].w2[2][0]	=	6.11381E+01	;
    WV[3].w2[3][0]	=	-6.89613E-02	;
    WV[3].w2[4][0]	=	1.54260E+01	;
    WV[3].w2[5][0]	=	-2.34379E+00	;
    WV[3].w2[6][0]	=	4.51307E+00	;
    WV[3].w2[7][0]	=	-3.88168E+01	;
    WV[3].w2[8][0]	=	-2.09305E-01	;
    WV[3].w2[9][0]	=	-2.28957E+01	;
    WV[3].w2[10][0]	=	1.83818E+01	;
    WV[3].w2[11][0]	=	-1.20588E+00	;
    WV[3].w2[12][0]	=	-1.07955E-01	;
    WV[3].w2[13][0]	=	-5.73910E-01	;
    WV[3].w2[14][0]	=	-5.30859E+00	;
    WV[3].w2[15][0]	=	1.70340E+01	; 
    nu_star_MAX[3] = 1.000000e+010;
    nu_star_MIN[3] = 1.000000e-015;
    D_star_MAX[3]  = 1.000000e+010;
    D_star_MIN[3]  = 1.000000e-010;
    GGGMAX[3]      = 1.000000e-000;
#endif
    
/*////////////////   RAX 375 data ////////////////*/
WV[4].omega[0] = -2.427730e+000;
WV[4].omega[1] = 1.804221e+000;
WV[4].omega[2] = -1.029778e+000;
WV[4].omega[3] = 5.036674e+000;
WV[4].omega[4] = -6.899825e-001;
WV[4].omega[5] = -3.527315e+000;
WV[4].omega[6] = -3.026513e+000;
WV[4].omega[7] = 1.555806e+000;
WV[4].omega[8] = -3.109531e+000;
WV[4].omega[9] = 1.691021e-001;
WV[4].omega[10] = -1.385853e+000;
WV[4].omega[11] = -2.976611e+000;
WV[4].omega[12] = 2.246877e+000;
WV[4].omega[13] = -2.080155e+000;
WV[4].omega[14] = 1.677169e+000;
WV[4].omega[15] = 2.755228e+000;
WV[4].omega[16] = 2.429726e-001;
WV[4].omega[17] = 1.588930e+000;
WV[4].omega[18] = -3.051072e+000;
WV[4].omega[19] = -9.666364e-001;
WV[4].omega[20] = -1.763031e+001;
WV[4].omega[21] = 1.017488e-001;
WV[4].omega[22] = -1.969543e+000;
WV[4].omega[23] = 2.106892e+000;
WV[4].omega[24] = -6.954959e-001;
WV[4].omega[25] = 4.075911e+000;
WV[4].omega[26] = 2.166805e+000;
WV[4].omega[27] = 6.305793e+000;
WV[4].omega[28] = 2.417978e+000;
WV[4].omega[29] = 2.941326e-001;
WV[4].omega[30] = -2.449323e+000;
WV[4].omega[31] = -1.759109e-001;
WV[4].omega[32] = -1.394707e+000;
WV[4].omega[33] = 2.742720e-001;
WV[4].omega[34] = -3.158252e-002;
WV[4].omega[35] = 6.296239e-001;
WV[4].omega[36] = 3.797693e-001;
WV[4].omega[37] = 1.808685e+000;
WV[4].omega[38] = 9.368261e-001;
WV[4].omega[39] = 2.009979e-001;
WV[4].omega[40] = 1.757464e+000;
WV[4].omega[41] = 9.626645e-001;
WV[4].omega[42] = -4.127616e-001;
WV[4].omega[43] = -2.198054e+000;
WV[4].omega[44] = -1.786763e-001;
WV[4].omega[45] = 1.418030e-002;
WV[4].omega[46] = 1.201606e-001;
WV[4].omega[47] = -2.329475e-001;
WV[4].omega[48] = 6.937171e-001;
WV[4].omega[49] = 6.877441e-001;
WV[4].omega[50] = -1.579221e+001;
WV[4].omega[51] = -1.387138e+000;
WV[4].omega[52] = 1.611863e-001;
WV[4].omega[53] = -1.167472e-001;
WV[4].omega[54] = 1.567574e-001;
WV[4].omega[55] = 1.849275e+000;
WV[4].omega[56] = -1.662574e-001;
WV[4].omega[57] = 7.713246e+000;
WV[4].omega[58] = -6.459472e-002;
WV[4].omega[59] = 1.745666e-001;
WV[4].omega[60] = -5.612760e+000;
WV[4].omega[61] = -1.316195e+001;
WV[4].omega[62] = 5.688894e+000;
WV[4].omega[63] = 1.061092e-001;
WV[4].omega[64] = -1.492889e+000;
WV[4].omega[65] = -3.445807e-001;
WV[4].omega[66] = 8.226556e-001;
WV[4].omega[67] = 1.234583e+001;
WV[4].omega[68] = 2.189386e+000;
WV[4].omega[69] = 3.335696e+000;
WV[4].omega[70] = -4.348153e-002;
WV[4].omega[71] = -2.258862e+000;
WV[4].omega[72] = -3.830243e+000;
WV[4].omega[73] = 1.303094e+001;
WV[4].omega[74] = 1.412948e+001;
WV[4].omega[75] = 3.834663e+000;
WV[4].w1[0][0] = -2.427730e+000;
WV[4].w1[0][1] = 1.804221e+000;
WV[4].w1[0][2] = -1.029778e+000;
WV[4].w1[0][3] = 5.036674e+000;
WV[4].w1[0][4] = -6.899825e-001;
WV[4].w1[0][5] = -3.527315e+000;
WV[4].w1[0][6] = -3.026513e+000;
WV[4].w1[0][7] = 1.555806e+000;
WV[4].w1[0][8] = -3.109531e+000;
WV[4].w1[0][9] = 1.691021e-001;
WV[4].w1[0][10] = -1.385853e+000;
WV[4].w1[0][11] = -2.976611e+000;
WV[4].w1[0][12] = 2.246877e+000;
WV[4].w1[0][13] = -2.080155e+000;
WV[4].w1[0][14] = 1.677169e+000;
WV[4].w1[1][0] = 2.755228e+000;
WV[4].w1[1][1] = 2.429726e-001;
WV[4].w1[1][2] = 1.588930e+000;
WV[4].w1[1][3] = -3.051072e+000;
WV[4].w1[1][4] = -9.666364e-001;
WV[4].w1[1][5] = -1.763031e+001;
WV[4].w1[1][6] = 1.017488e-001;
WV[4].w1[1][7] = -1.969543e+000;
WV[4].w1[1][8] = 2.106892e+000;
WV[4].w1[1][9] = -6.954959e-001;
WV[4].w1[1][10] = 4.075911e+000;
WV[4].w1[1][11] = 2.166805e+000;
WV[4].w1[1][12] = 6.305793e+000;
WV[4].w1[1][13] = 2.417978e+000;
WV[4].w1[1][14] = 2.941326e-001;
WV[4].w1[2][0] = -2.449323e+000;
WV[4].w1[2][1] = -1.759109e-001;
WV[4].w1[2][2] = -1.394707e+000;
WV[4].w1[2][3] = 2.742720e-001;
WV[4].w1[2][4] = -3.158252e-002;
WV[4].w1[2][5] = 6.296239e-001;
WV[4].w1[2][6] = 3.797693e-001;
WV[4].w1[2][7] = 1.808685e+000;
WV[4].w1[2][8] = 9.368261e-001;
WV[4].w1[2][9] = 2.009979e-001;
WV[4].w1[2][10] = 1.757464e+000;
WV[4].w1[2][11] = 9.626645e-001;
WV[4].w1[2][12] = -4.127616e-001;
WV[4].w1[2][13] = -2.198054e+000;
WV[4].w1[2][14] = -1.786763e-001;
WV[4].w1[3][0] = 1.418030e-002;
WV[4].w1[3][1] = 1.201606e-001;
WV[4].w1[3][2] = -2.329475e-001;
WV[4].w1[3][3] = 6.937171e-001;
WV[4].w1[3][4] = 6.877441e-001;
WV[4].w1[3][5] = -1.579221e+001;
WV[4].w1[3][6] = -1.387138e+000;
WV[4].w1[3][7] = 1.611863e-001;
WV[4].w1[3][8] = -1.167472e-001;
WV[4].w1[3][9] = 1.567574e-001;
WV[4].w1[3][10] = 1.849275e+000;
WV[4].w1[3][11] = -1.662574e-001;
WV[4].w1[3][12] = 7.713246e+000;
WV[4].w1[3][13] = -6.459472e-002;
WV[4].w1[3][14] = 1.745666e-001;
WV[4].w2[0][0] = -5.612760e+000;
WV[4].w2[1][0] = -1.316195e+001;
WV[4].w2[2][0] = 5.688894e+000;
WV[4].w2[3][0] = 1.061092e-001;
WV[4].w2[4][0] = -1.492889e+000;
WV[4].w2[5][0] = -3.445807e-001;
WV[4].w2[6][0] = 8.226556e-001;
WV[4].w2[7][0] = 1.234583e+001;
WV[4].w2[8][0] = 2.189386e+000;
WV[4].w2[9][0] = 3.335696e+000;
WV[4].w2[10][0] = -4.348153e-002;
WV[4].w2[11][0] = -2.258862e+000;
WV[4].w2[12][0] = -3.830243e+000;
WV[4].w2[13][0] = 1.303094e+001;
WV[4].w2[14][0] = 1.412948e+001;
WV[4].w2[15][0] = 3.834663e+000;
nu_star_MAX[4] = 1.000000e+003;
nu_star_MIN[4] = 1.000000e-006;
D_star_MAX[4]  = 1.000000e+004;
D_star_MIN[4]  = 1.000000e-004;
GGGMAX[4]      = 3.000000e-001;

/*////////////////   RAX 390 data ////////////////*/
WV[5].omega[0] = 	-8.807072E-01	;
WV[5].omega[1] = 	8.007746E-01	;
WV[5].omega[2] = 	4.147654E+00	;
WV[5].omega[3] = 	-1.931057E+00	;
WV[5].omega[4] = 	-3.646674E+00	;
WV[5].omega[5] = 	-2.923838E+00	;
WV[5].omega[6] = 	-2.191935E+00	;
WV[5].omega[7] = 	2.689461E+00	;
WV[5].omega[8] = 	-1.715148E+00	;
WV[5].omega[9] = 	7.871025E+00	;
WV[5].omega[10] =	1.204518E+00	;
WV[5].omega[11] =	3.705568E+00	;
WV[5].omega[12] =	-1.276198E+00	;
WV[5].omega[13] =	-2.129330E+00	;
WV[5].omega[14] =	1.124981E+00	;
WV[5].omega[15] =	-8.125661E-01	;
WV[5].omega[16] =	-1.891950E+00	;
WV[5].omega[17] =	5.010538E-01	;
WV[5].omega[18] =	-2.150450E+00	;
WV[5].omega[19] =	6.181748E+00	;
WV[5].omega[20] =	3.029499E+00	;
WV[5].omega[21] =	1.523311E+00	;
WV[5].omega[22] =	-2.707538E+00	;
WV[5].omega[23] =	-4.919353E-01	;
WV[5].omega[24] =	6.100653E+00	;
WV[5].omega[25] =	-7.257975E-01	;
WV[5].omega[26] =	-6.445360E+00	;
WV[5].omega[27] =	5.474779E-01	;
WV[5].omega[28] =	1.497906E+00	;
WV[5].omega[29] =	1.241128E+00	;
WV[5].omega[30] =	-5.467468E-02	;
WV[5].omega[31] =	6.346076E-02	;
WV[5].omega[32] =	-2.222877E-01	;
WV[5].omega[33] =	2.416387E+00	;
WV[5].omega[34] =	5.664591E-01	;
WV[5].omega[35] =	4.633641E-01	;
WV[5].omega[36] =	4.417510E+00	;
WV[5].omega[37] =	-4.280496E-01	;
WV[5].omega[38] =	-5.743868E-02	;
WV[5].omega[39] =	-9.649082E-01	;
WV[5].omega[40] =	1.174685E-01	;
WV[5].omega[41] =	-5.742852E-01	;
WV[5].omega[42] =	-5.477763E-01	;
WV[5].omega[43] =	3.884035E+00	;
WV[5].omega[44] =	-3.398475E-01	;
WV[5].omega[45] =	9.127401E-01	;
WV[5].omega[46] =	1.726610E+00	;
WV[5].omega[47] =	-1.275442E+00	;
WV[5].omega[48] =	-4.225579E+00	;
WV[5].omega[49] =	1.092599E+00	;
WV[5].omega[50] =	1.283578E+00	;
WV[5].omega[51] =	-8.583820E-01	;
WV[5].omega[52] =	-1.050956E+00	;
WV[5].omega[53] =	1.524835E+00	;
WV[5].omega[54] =	1.273782E+01	;
WV[5].omega[55] =	1.213281E+00	;
WV[5].omega[56] =	-1.148391E+00	;
WV[5].omega[57] =	-1.883974E+00	;
WV[5].omega[58] =	-6.396538E-01	;
WV[5].omega[59] =	2.714039E+00	;
WV[5].omega[60] =	1.206796E+00	;
WV[5].omega[61] =	-2.545047E+00	;
WV[5].omega[62] =	1.464093E-01	;
WV[5].omega[63] =	-1.014409E-01	;
WV[5].omega[64] =	-3.145115E+00	;
WV[5].omega[65] =	-2.097942E+00	;
WV[5].omega[66] =	3.688845E-01	;
WV[5].omega[67] =	-2.762660E+00	;
WV[5].omega[68] =	-1.575968E+00	;
WV[5].omega[69] =	-2.136212E+00	;
WV[5].omega[70] =	3.143207E+00	;
WV[5].omega[71] =	-2.817999E+00	;
WV[5].omega[72] =	1.636712E+00	;
WV[5].omega[73] =	-4.067499E-01	;
WV[5].omega[74] =	-1.963729E+00	;
WV[5].omega[75] =	5.751423E+00	;
WV[5].w1[0][0] = 	-8.807072E-01	;
WV[5].w1[0][1] = 	8.007746E-01	;
WV[5].w1[0][2] = 	4.147654E+00	;
WV[5].w1[0][3] = 	-1.931057E+00	;
WV[5].w1[0][4] = 	-3.646674E+00	;
WV[5].w1[0][5] = 	-2.923838E+00	;
WV[5].w1[0][6] = 	-2.191935E+00	;
WV[5].w1[0][7] = 	2.689461E+00	;
WV[5].w1[0][8] = 	-1.715148E+00	;
WV[5].w1[0][9] = 	7.871025E+00	;
WV[5].w1[0][10] =	1.204518E+00	;
WV[5].w1[0][11] =	3.705568E+00	;
WV[5].w1[0][12] =	-1.276198E+00	;
WV[5].w1[0][13] =	-2.129330E+00	;
WV[5].w1[0][14] =	1.124981E+00	;
WV[5].w1[1][0] = 	-8.125661E-01	;
WV[5].w1[1][1] = 	-1.891950E+00	;
WV[5].w1[1][2] = 	5.010538E-01	;
WV[5].w1[1][3] = 	-2.150450E+00	;
WV[5].w1[1][4] = 	6.181748E+00	;
WV[5].w1[1][5] = 	3.029499E+00	;
WV[5].w1[1][6] = 	1.523311E+00	;
WV[5].w1[1][7] = 	-2.707538E+00	;
WV[5].w1[1][8] = 	-4.919353E-01	;
WV[5].w1[1][9] = 	6.100653E+00	;
WV[5].w1[1][10] =	-7.257975E-01	;
WV[5].w1[1][11] =	-6.445360E+00	;
WV[5].w1[1][12] =	5.474779E-01	;
WV[5].w1[1][13] =	1.497906E+00	;
WV[5].w1[1][14] =	1.241128E+00	;
WV[5].w1[2][0] = 	-5.467468E-02	;
WV[5].w1[2][1] = 	6.346076E-02	;
WV[5].w1[2][2] = 	-2.222877E-01	;
WV[5].w1[2][3] = 	2.416387E+00	;
WV[5].w1[2][4] = 	5.664591E-01	;
WV[5].w1[2][5] = 	4.633641E-01	;
WV[5].w1[2][6] = 	4.417510E+00	;
WV[5].w1[2][7] = 	-4.280496E-01	;
WV[5].w1[2][8] = 	-5.743868E-02	;
WV[5].w1[2][9] = 	-9.649082E-01	;
WV[5].w1[2][10] =	1.174685E-01	;
WV[5].w1[2][11] =	-5.742852E-01	;
WV[5].w1[2][12] =	-5.477763E-01	;
WV[5].w1[2][13] =	3.884035E+00	;
WV[5].w1[2][14] =	-3.398475E-01	;
WV[5].w1[3][0] = 	9.127401E-01	;
WV[5].w1[3][1] = 	1.726610E+00	;
WV[5].w1[3][2] = 	-1.275442E+00	;
WV[5].w1[3][3] = 	-4.225579E+00	;
WV[5].w1[3][4] = 	1.092599E+00	;
WV[5].w1[3][5] = 	1.283578E+00	;
WV[5].w1[3][6] = 	-8.583820E-01	;
WV[5].w1[3][7] = 	-1.050956E+00	;
WV[5].w1[3][8] = 	1.524835E+00	;
WV[5].w1[3][9] = 	1.273782E+01	;
WV[5].w1[3][10] =	1.213281E+00	;
WV[5].w1[3][11] =	-1.148391E+00	;
WV[5].w1[3][12] =	-1.883974E+00	;
WV[5].w1[3][13] =	-6.396538E-01	;
WV[5].w1[3][14] =	2.714039E+00	;
WV[5].w2[0][0] = 	1.206796E+00	;
WV[5].w2[1][0] = 	-2.545047E+00	;
WV[5].w2[2][0] = 	1.464093E-01	;
WV[5].w2[3][0] = 	-1.014409E-01	;
WV[5].w2[4][0] = 	-3.145115E+00	;
WV[5].w2[5][0] = 	-2.097942E+00	;
WV[5].w2[6][0] = 	3.688845E-01	;
WV[5].w2[7][0] = 	-2.762660E+00	;
WV[5].w2[8][0] = 	-1.575968E+00	;
WV[5].w2[9][0] = 	-2.136212E+00	;
WV[5].w2[10][0] =	3.143207E+00	;
WV[5].w2[11][0] =	-2.817999E+00	;
WV[5].w2[12][0] =	1.636712E+00	;
WV[5].w2[13][0] =	-4.067499E-01	;
WV[5].w2[14][0] =	-1.963729E+00	;
WV[5].w2[15][0] =	5.751423E+00	;
nu_star_MAX[5] = 1.000000e+003;
nu_star_MIN[5] = 1.000000e-006;
D_star_MAX[5]  = 1.000000e+005;
D_star_MIN[5]  = 1.000000e-004;
GGGMAX[5]      = 3.120000e-001;

/*////////////////   RAX 360 beta 1% data ////////////////*/		
/* 2009.09.30 このNNWをシヨウした際に、計算結果においてD1[1]が10^5のオーダになっている. */
/* 要確認*/
//WV[6].omega[0] = 	-2.213E+00	;
//WV[6].omega[1] = 	2.594E+00	;
//WV[6].omega[2] = 	6.015E+00	;
//WV[6].omega[3] = 	-2.188E+00	;
//WV[6].omega[4] = 	-2.445E+00	;
//WV[6].omega[5] = 	1.153E+00	;
//WV[6].omega[6] = 	3.112E+00	;
//WV[6].omega[7] = 	-4.991E+00	;
//WV[6].omega[8] = 	-5.472E+00	;
//WV[6].omega[9] = 	-3.456E+00	;
//WV[6].omega[10] =	-4.467E+00	;
//WV[6].omega[11] =	1.343E+00	;
//WV[6].omega[12] =	8.026E-01	;
//WV[6].omega[13] =	-1.565E+00	;
//WV[6].omega[14] =	3.084E+00	;
//WV[6].omega[15] =	-9.565E+00	;
//WV[6].omega[16] =	-3.467E+00	;
//WV[6].omega[17] =	1.070E+00	;
//WV[6].omega[18] =	1.575E+00	;
//WV[6].omega[19] =	1.511E+00	;
//WV[6].omega[20] =	-3.242E+01	;
//WV[6].omega[21] =	-2.924E+00	;
//WV[6].omega[22] =	-3.472E-01	;
//WV[6].omega[23] =	3.558E+00	;
//WV[6].omega[24] =	2.720E+00	;
//WV[6].omega[25] =	-9.133E+00	;
//WV[6].omega[26] =	-3.066E+01	;
//WV[6].omega[27] =	-3.569E+01	;
//WV[6].omega[28] =	1.802E+00	;
//WV[6].omega[29] =	-3.153E+00	;
//WV[6].omega[30] =	-1.170E+00	;
//WV[6].omega[31] =	7.624E-01	;
//WV[6].omega[32] =	8.930E-01	;
//WV[6].omega[33] =	-4.413E-01	;
//WV[6].omega[34] =	-4.523E-01	;
//WV[6].omega[35] =	4.535E+00	;
//WV[6].omega[36] =	-9.128E-01	;
//WV[6].omega[37] =	7.642E-02	;
//WV[6].omega[38] =	-1.158E+00	;
//WV[6].omega[39] =	8.256E-01	;
//WV[6].omega[40] =	-1.914E-01	;
//WV[6].omega[41] =	4.282E+00	;
//WV[6].omega[42] =	5.009E+00	;
//WV[6].omega[43] =	-4.300E-01	;
//WV[6].omega[44] =	7.624E-01	;
//WV[6].omega[45] =	-1.018E+01	;
//WV[6].omega[46] =	-4.088E-02	;
//WV[6].omega[47] =	2.826E+00	;
//WV[6].omega[48] =	6.650E-01	;
//WV[6].omega[49] =	7.913E-01	;
//WV[6].omega[50] =	-7.853E+00	;
//WV[6].omega[51] =	-1.443E+00	;
//WV[6].omega[52] =	2.758E+00	;
//WV[6].omega[53] =	-7.179E-01	;
//WV[6].omega[54] =	1.547E+00	;
//WV[6].omega[55] =	-9.954E+00	;
//WV[6].omega[56] =	-7.498E+00	;
//WV[6].omega[57] =	-8.513E+00	;
//WV[6].omega[58] =	-5.493E-01	;
//WV[6].omega[59] =	1.361E+00	;
//WV[6].omega[60] =	-1.049E+00	;
//WV[6].omega[61] =	2.483E+00	;
//WV[6].omega[62] =	-3.106E-01	;
//WV[6].omega[63] =	-6.706E+00	;
//WV[6].omega[64] =	4.477E+00	;
//WV[6].omega[65] =	-1.293E+01	;
//WV[6].omega[66] =	5.264E-01	;
//WV[6].omega[67] =	-1.600E-01	;
//WV[6].omega[68] =	-3.956E-01	;
//WV[6].omega[69] =	5.861E-01	;
//WV[6].omega[70] =	1.028E+00	;
//WV[6].omega[71] =	8.120E+00	;
//WV[6].omega[72] =	4.816E+00	;
//WV[6].omega[73] =	1.280E+01	;
//WV[6].omega[74] =	6.018E+00	;
//WV[6].omega[75] =	2.188E+00	;
//WV[6].w1[0][0] = 	-2.213E+00	;
//WV[6].w1[0][1] = 	2.594E+00	;
//WV[6].w1[0][2] = 	6.015E+00	;
//WV[6].w1[0][3] = 	-2.188E+00	;
//WV[6].w1[0][4] = 	-2.445E+00	;
//WV[6].w1[0][5] = 	1.153E+00	;
//WV[6].w1[0][6] = 	3.112E+00	;
//WV[6].w1[0][7] = 	-4.991E+00	;
//WV[6].w1[0][8] = 	-5.472E+00	;
//WV[6].w1[0][9] = 	-3.456E+00	;
//WV[6].w1[0][10] =	-4.467E+00	;
//WV[6].w1[0][11] =	1.343E+00	;
//WV[6].w1[0][12] =	8.026E-01	;
//WV[6].w1[0][13] =	-1.565E+00	;
//WV[6].w1[0][14] =	3.084E+00	;
//WV[6].w1[1][0] = 	-9.565E+00	;
//WV[6].w1[1][1] = 	-3.467E+00	;
//WV[6].w1[1][2] = 	1.070E+00	;
//WV[6].w1[1][3] = 	1.575E+00	;
//WV[6].w1[1][4] = 	1.511E+00	;
//WV[6].w1[1][5] = 	-3.242E+01	;
//WV[6].w1[1][6] = 	-2.924E+00	;
//WV[6].w1[1][7] = 	-3.472E-01	;
//WV[6].w1[1][8] = 	3.558E+00	;
//WV[6].w1[1][9] = 	2.720E+00	;
//WV[6].w1[1][10] =	-9.133E+00	;
//WV[6].w1[1][11] =	-3.066E+01	;
//WV[6].w1[1][12] =	-3.569E+01	;
//WV[6].w1[1][13] =	1.802E+00	;
//WV[6].w1[1][14] =	-3.153E+00	;
//WV[6].w1[2][0] = 	-1.170E+00	;
//WV[6].w1[2][1] = 	7.624E-01	;
//WV[6].w1[2][2] = 	8.930E-01	;
//WV[6].w1[2][3] = 	-4.413E-01	;
//WV[6].w1[2][4] = 	-4.523E-01	;
//WV[6].w1[2][5] = 	4.535E+00	;
//WV[6].w1[2][6] = 	-9.128E-01	;
//WV[6].w1[2][7] = 	7.642E-02	;
//WV[6].w1[2][8] = 	-1.158E+00	;
//WV[6].w1[2][9] = 	8.256E-01	;
//WV[6].w1[2][10] =	-1.914E-01	;
//WV[6].w1[2][11] =	4.282E+00	;
//WV[6].w1[2][12] =	5.009E+00	;
//WV[6].w1[2][13] =	-4.300E-01	;
//WV[6].w1[2][14] =	7.624E-01	;
//WV[6].w1[3][0] = 	-1.018E+01	;
//WV[6].w1[3][1] = 	-4.088E-02	;
//WV[6].w1[3][2] = 	2.826E+00	;
//WV[6].w1[3][3] = 	6.650E-01	;
//WV[6].w1[3][4] = 	7.913E-01	;
//WV[6].w1[3][5] = 	-7.853E+00	;
//WV[6].w1[3][6] = 	-1.443E+00	;
//WV[6].w1[3][7] = 	2.758E+00	;
//WV[6].w1[3][8] = 	-7.179E-01	;
//WV[6].w1[3][9] = 	1.547E+00	;
//WV[6].w1[3][10] =	-9.954E+00	;
//WV[6].w1[3][11] =	-7.498E+00	;
//WV[6].w1[3][12] =	-8.513E+00	;
//WV[6].w1[3][13] =	-5.493E-01	;
//WV[6].w1[3][14] =	1.361E+00	;
//WV[6].w2[0][0] = 	-1.049E+00	;
//WV[6].w2[1][0] = 	2.483E+00	;
//WV[6].w2[2][0] = 	-3.106E-01	;
//WV[6].w2[3][0] = 	-6.706E+00	;
//WV[6].w2[4][0] = 	4.477E+00	;
//WV[6].w2[5][0] = 	-1.293E+01	;
//WV[6].w2[6][0] = 	5.264E-01	;
//WV[6].w2[7][0] = 	-1.600E-01	;
//WV[6].w2[8][0] = 	-3.956E-01	;
//WV[6].w2[9][0] = 	5.861E-01	;
//WV[6].w2[10][0] =	1.028E+00	;
//WV[6].w2[11][0] =	8.120E+00	;
//WV[6].w2[12][0] =	4.816E+00	;
//WV[6].w2[13][0] =	1.280E+01	;
//WV[6].w2[14][0] =	6.018E+00	;
//WV[6].w2[15][0] =	2.188E+00	;
//nu_star_MAX[6] = 	1.000E+03	;
//nu_star_MIN[6] = 	1.000E-10	;
//D_star_MAX[6]  = 	1.000E+08	;
//D_star_MIN[6]  = 	1.000E-08	;
//GGGMAX[6]      = 	3.160E-01	;


//	/*   RAX 360 beta 1% data */		
//	/* 2009.10.02 2009.7月に更新版（？）と思われるNNW */		
//	/* 要確認*/		
//	WV[6].omega[0] = 	3.90E+00	;
//	WV[6].omega[1] = 	-8.64E+00	;
//	WV[6].omega[2] = 	2.54E+00	;
//	WV[6].omega[3] = 	2.70E+00	;
//	WV[6].omega[4] = 	-3.34E+00	;
//	WV[6].omega[5] = 	3.13E+00	;
//	WV[6].omega[6] = 	2.09E+02	;
//	WV[6].omega[7] = 	-1.30E+00	;
//	WV[6].omega[8] = 	5.71E-01	;
//	WV[6].omega[9] = 	-7.24E-02	;
//	WV[6].omega[10] =	1.32E-01	;
//	WV[6].omega[11] =	-7.02E+00	;
//	WV[6].omega[12] =	-4.81E+00	;
//	WV[6].omega[13] =	-3.85E+01	;
//	WV[6].omega[14] =	1.67E+00	;
//	WV[6].omega[15] =	3.87E+01	;
//	WV[6].omega[16] =	9.68E-02	;
//	WV[6].omega[17] =	-9.97E-01	;
//	WV[6].omega[18] =	-1.07E+00	;
//	WV[6].omega[19] =	-5.68E-01	;
//	WV[6].omega[20] =	-1.12E+00	;
//	WV[6].omega[21] =	-8.34E+00	;
//	WV[6].omega[22] =	-1.26E+00	;
//	WV[6].omega[23] =	1.15E+00	;
//	WV[6].omega[24] =	-3.09E+00	;
//	WV[6].omega[25] =	3.07E+00	;
//	WV[6].omega[26] =	-1.28E+02	;
//	WV[6].omega[27] =	2.72E+00	;
//	WV[6].omega[28] =	2.42E+01	;
//	WV[6].omega[29] =	1.35E+00	;
//	WV[6].omega[30] =	2.48E+00	;
//	WV[6].omega[31] =	-5.99E-03	;
//	WV[6].omega[32] =	-1.05E-01	;
//	WV[6].omega[33] =	-9.40E-02	;
//	WV[6].omega[34] =	6.02E-01	;
//	WV[6].omega[35] =	-1.88E-01	;
//	WV[6].omega[36] =	-2.43E+01	;
//	WV[6].omega[37] =	3.43E-01	;
//	WV[6].omega[38] =	-2.64E-01	;
//	WV[6].omega[39] =	-6.78E-01	;
//	WV[6].omega[40] =	6.57E-01	;
//	WV[6].omega[41] =	7.96E-01	;
//	WV[6].omega[42] =	8.34E-01	;
//	WV[6].omega[43] =	2.77E+00	;
//	WV[6].omega[44] =	-3.80E-01	;
//	WV[6].omega[45] =	3.39E+01	;
//	WV[6].omega[46] =	3.84E+00	;
//	WV[6].omega[47] =	-3.06E-01	;
//	WV[6].omega[48] =	-2.73E-01	;
//	WV[6].omega[49] =	-1.78E+00	;
//	WV[6].omega[50] =	-1.05E+00	;
//	WV[6].omega[51] =	9.24E+01	;
//	WV[6].omega[52] =	1.73E-01	;
//	WV[6].omega[53] =	8.61E-02	;
//	WV[6].omega[54] =	-9.89E-01	;
//	WV[6].omega[55] =	9.56E-01	;
//	WV[6].omega[56] =	-1.05E+02	;
//	WV[6].omega[57] =	6.04E-01	;
//	WV[6].omega[58] =	-7.26E+00	;
//	WV[6].omega[59] =	-2.88E-01	;
//	WV[6].omega[60] =	1.05E-01	;
//	WV[6].omega[61] =	-1.22E-01	;
//	WV[6].omega[62] =	-1.93E+01	;
//	WV[6].omega[63] =	1.60E+01	;
//	WV[6].omega[64] =	-7.22E-01	;
//	WV[6].omega[65] =	1.78E+00	;
//	WV[6].omega[66] =	-3.52E-02	;
//	WV[6].omega[67] =	-9.66E+00	;
//	WV[6].omega[68] =	-4.80E+00	;
//	WV[6].omega[69] =	2.13E+00	;
//	WV[6].omega[70] =	2.17E+00	;
//	WV[6].omega[71] =	3.64E-01	;
//	WV[6].omega[72] =	-3.01E-01	;
//	WV[6].omega[73] =	-4.53E-02	;
//	WV[6].omega[74] =	-5.33E+00	;
//	WV[6].omega[75] =	2.83E-01	;
//	WV[6].w1[0][0] = 	3.90E+00	;
//	WV[6].w1[0][1] = 	-8.64E+00	;
//	WV[6].w1[0][2] = 	2.54E+00	;
//	WV[6].w1[0][3] = 	2.70E+00	;
//	WV[6].w1[0][4] = 	-3.34E+00	;
//	WV[6].w1[0][5] = 	3.13E+00	;
//	WV[6].w1[0][6] = 	2.09E+02	;
//	WV[6].w1[0][7] = 	-1.30E+00	;
//	WV[6].w1[0][8] = 	5.71E-01	;
//	WV[6].w1[0][9] = 	-7.24E-02	;
//	WV[6].w1[0][10] =	1.32E-01	;
//	WV[6].w1[0][11] =	-7.02E+00	;
//	WV[6].w1[0][12] =	-4.81E+00	;
//	WV[6].w1[0][13] =	-3.85E+01	;
//	WV[6].w1[0][14] =	1.67E+00	;
//	WV[6].w1[1][0] = 	3.87E+01	;
//	WV[6].w1[1][1] = 	9.68E-02	;
//	WV[6].w1[1][2] = 	-9.97E-01	;
//	WV[6].w1[1][3] = 	-1.07E+00	;
//	WV[6].w1[1][4] = 	-5.68E-01	;
//	WV[6].w1[1][5] = 	-1.12E+00	;
//	WV[6].w1[1][6] = 	-8.34E+00	;
//	WV[6].w1[1][7] = 	-1.26E+00	;
//	WV[6].w1[1][8] = 	1.15E+00	;
//	WV[6].w1[1][9] = 	-3.09E+00	;
//	WV[6].w1[1][10] =	3.07E+00	;
//	WV[6].w1[1][11] =	-1.28E+02	;
//	WV[6].w1[1][12] =	2.72E+00	;
//	WV[6].w1[1][13] =	2.42E+01	;
//	WV[6].w1[1][14] =	1.35E+00	;
//	WV[6].w1[2][0] = 	2.48E+00	;
//	WV[6].w1[2][1] = 	-5.99E-03	;
//	WV[6].w1[2][2] = 	-1.05E-01	;
//	WV[6].w1[2][3] = 	-9.40E-02	;
//	WV[6].w1[2][4] = 	6.02E-01	;
//	WV[6].w1[2][5] = 	-1.88E-01	;
//	WV[6].w1[2][6] = 	-2.43E+01	;
//	WV[6].w1[2][7] = 	3.43E-01	;
//	WV[6].w1[2][8] = 	-2.64E-01	;
//	WV[6].w1[2][9] = 	-6.78E-01	;
//	WV[6].w1[2][10] =	6.57E-01	;
//	WV[6].w1[2][11] =	7.96E-01	;
//	WV[6].w1[2][12] =	8.34E-01	;
//	WV[6].w1[2][13] =	2.77E+00	;
//	WV[6].w1[2][14] =	-3.80E-01	;
//	WV[6].w1[3][0] = 	3.39E+01	;
//	WV[6].w1[3][1] = 	3.84E+00	;
//	WV[6].w1[3][2] = 	-3.06E-01	;
//	WV[6].w1[3][3] = 	-2.73E-01	;
//	WV[6].w1[3][4] = 	-1.78E+00	;
//	WV[6].w1[3][5] = 	-1.05E+00	;
//	WV[6].w1[3][6] = 	9.24E+01	;
//	WV[6].w1[3][7] = 	1.73E-01	;
//	WV[6].w1[3][8] = 	8.61E-02	;
//	WV[6].w1[3][9] = 	-9.89E-01	;
//	WV[6].w1[3][10] =	9.56E-01	;
//	WV[6].w1[3][11] =	-1.05E+02	;
//	WV[6].w1[3][12] =	6.04E-01	;
//	WV[6].w1[3][13] =	-7.26E+00	;
//	WV[6].w1[3][14] =	-2.88E-01	;
//	WV[6].w2[0][0] = 	1.05E-01	;
//	WV[6].w2[1][0] = 	-1.22E-01	;
//	WV[6].w2[2][0] = 	-1.93E+01	;
//	WV[6].w2[3][0] = 	1.60E+01	;
//	WV[6].w2[4][0] = 	-7.22E-01	;
//	WV[6].w2[5][0] = 	1.78E+00	;
//	WV[6].w2[6][0] = 	-3.52E-02	;
//	WV[6].w2[7][0] = 	-9.66E+00	;
//	WV[6].w2[8][0] = 	-4.80E+00	;
//	WV[6].w2[9][0] = 	2.13E+00	;
//	WV[6].w2[10][0] =	2.17E+00	;
//	WV[6].w2[11][0] =	3.64E-01	;
//	WV[6].w2[12][0] =	-3.01E-01	;
//	WV[6].w2[13][0] =	-4.53E-02	;
//	WV[6].w2[14][0] =	-5.33E+00	;
//	WV[6].w2[15][0] =	2.83E-01	;
//	nu_star_MAX[6] = 	1.00E+05	;
//	nu_star_MIN[6] = 	1.00E-12	;
//	D_star_MAX[6]  = 	1.00E+10	;
//	D_star_MIN[6]  = 	1.00E-10	;
//	GGGMAX[6]      = 	3.16E-01	;

/*   RAX 360 beta 1% data */		

//	WV[6].omega[0] = 	-5.84E+00	;
//	WV[6].omega[1] = 	-4.96E+00	;
//	WV[6].omega[2] = 	-1.36E+00	;
//	WV[6].omega[3] = 	3.13E+00	;
//	WV[6].omega[4] = 	4.61E+01	;
//	WV[6].omega[5] = 	1.02E+01	;
//	WV[6].omega[6] = 	3.31E+00	;
//	WV[6].omega[7] = 	-4.20E+00	;
//	WV[6].omega[8] = 	7.63E+01	;
//	WV[6].omega[9] = 	4.97E+00	;
//	WV[6].omega[10] =	7.58E+01	;
//	WV[6].omega[11] =	2.40E+00	;
//	WV[6].omega[12] =	1.13E+00	;
//	WV[6].omega[13] =	9.17E+01	;
//	WV[6].omega[14] =	2.31E+00	;
//	WV[6].omega[15] =	3.61E+00	;
//	WV[6].omega[16] =	-2.60E+02	;
//	WV[6].omega[17] =	1.79E-01	;
//	WV[6].omega[18] =	-6.22E-01	;
//	WV[6].omega[19] =	-2.65E+01	;
//	WV[6].omega[20] =	2.10E-01	;
//	WV[6].omega[21] =	-7.77E-01	;
//	WV[6].omega[22] =	6.58E-01	;
//	WV[6].omega[23] =	-3.34E+01	;
//	WV[6].omega[24] =	1.05E+02	;
//	WV[6].omega[25] =	-3.32E+01	;
//	WV[6].omega[26] =	1.37E+00	;
//	WV[6].omega[27] =	-1.89E-01	;
//	WV[6].omega[28] =	-4.90E+01	;
//	WV[6].omega[29] =	1.38E+00	;
//	WV[6].omega[30] =	1.43E+00	;
//	WV[6].omega[31] =	9.81E-01	;
//	WV[6].omega[32] =	2.22E-01	;
//	WV[6].omega[33] =	-3.44E-01	;
//	WV[6].omega[34] =	9.45E+01	;
//	WV[6].omega[35] =	1.57E-01	;
//	WV[6].omega[36] =	-4.15E-01	;
//	WV[6].omega[37] =	6.53E-02	;
//	WV[6].omega[38] =	2.03E+02	;
//	WV[6].omega[39] =	-9.82E-01	;
//	WV[6].omega[40] =	1.92E+02	;
//	WV[6].omega[41] =	5.48E-01	;
//	WV[6].omega[42] =	-2.80E-01	;
//	WV[6].omega[43] =	-8.65E+00	;
//	WV[6].omega[44] =	3.60E-01	;
//	WV[6].omega[45] =	-6.39E-01	;
//	WV[6].omega[46] =	-2.10E+02	;
//	WV[6].omega[47] =	-1.01E+00	;
//	WV[6].omega[48] =	5.74E-01	;
//	WV[6].omega[49] =	3.72E+01	;
//	WV[6].omega[50] =	-4.36E+00	;
//	WV[6].omega[51] =	5.86E-01	;
//	WV[6].omega[52] =	-4.19E-01	;
//	WV[6].omega[53] =	2.71E+01	;
//	WV[6].omega[54] =	8.58E+01	;
//	WV[6].omega[55] =	2.34E+01	;
//	WV[6].omega[56] =	3.36E+00	;
//	WV[6].omega[57] =	1.44E+00	;
//	WV[6].omega[58] =	-3.52E+01	;
//	WV[6].omega[59] =	3.49E+00	;
//	WV[6].omega[60] =	2.62E-01	;
//	WV[6].omega[61] =	2.61E+01	;
//	WV[6].omega[62] =	-9.65E+00	;
//	WV[6].omega[63] =	-1.37E+01	;
//	WV[6].omega[64] =	-1.35E-03	;
//	WV[6].omega[65] =	7.35E-02	;
//	WV[6].omega[66] =	1.11E+01	;
//	WV[6].omega[67] =	-1.06E+00	;
//	WV[6].omega[68] =	-1.12E+00	;
//	WV[6].omega[69] =	2.58E+01	;
//	WV[6].omega[70] =	1.12E+00	;
//	WV[6].omega[71] =	8.60E+00	;
//	WV[6].omega[72] =	-1.09E+01	;
//	WV[6].omega[73] =	-4.49E-03	;
//	WV[6].omega[74] =	-1.46E+01	;
//	WV[6].omega[75] =	9.54E+00	;
//	WV[6].w1[0][0] = 	-5.84E+00	;
//	WV[6].w1[0][1] = 	-4.96E+00	;
//	WV[6].w1[0][2] = 	-1.36E+00	;
//	WV[6].w1[0][3] = 	3.13E+00	;
//	WV[6].w1[0][4] = 	4.61E+01	;
//	WV[6].w1[0][5] = 	1.02E+01	;
//	WV[6].w1[0][6] = 	3.31E+00	;
//	WV[6].w1[0][7] = 	-4.20E+00	;
//	WV[6].w1[0][8] = 	7.63E+01	;
//	WV[6].w1[0][9] = 	4.97E+00	;
//	WV[6].w1[0][10] =	7.58E+01	;
//	WV[6].w1[0][11] =	2.40E+00	;
//	WV[6].w1[0][12] =	1.13E+00	;
//	WV[6].w1[0][13] =	9.17E+01	;
//	WV[6].w1[0][14] =	2.31E+00	;
//	WV[6].w1[1][0] = 	3.61E+00	;
//	WV[6].w1[1][1] = 	-2.60E+02	;
//	WV[6].w1[1][2] = 	1.79E-01	;
//	WV[6].w1[1][3] = 	-6.22E-01	;
//	WV[6].w1[1][4] = 	-2.65E+01	;
//	WV[6].w1[1][5] = 	2.10E-01	;
//	WV[6].w1[1][6] = 	-7.77E-01	;
//	WV[6].w1[1][7] = 	6.58E-01	;
//	WV[6].w1[1][8] = 	-3.34E+01	;
//	WV[6].w1[1][9] = 	1.05E+02	;
//	WV[6].w1[1][10] =	-3.32E+01	;
//	WV[6].w1[1][11] =	1.37E+00	;
//	WV[6].w1[1][12] =	-1.89E-01	;
//	WV[6].w1[1][13] =	-4.90E+01	;
//	WV[6].w1[1][14] =	1.38E+00	;
//	WV[6].w1[2][0] = 	1.43E+00	;
//	WV[6].w1[2][1] = 	9.81E-01	;
//	WV[6].w1[2][2] = 	2.22E-01	;
//	WV[6].w1[2][3] = 	-3.44E-01	;
//	WV[6].w1[2][4] = 	9.45E+01	;
//	WV[6].w1[2][5] = 	1.57E-01	;
//	WV[6].w1[2][6] = 	-4.15E-01	;
//	WV[6].w1[2][7] = 	6.53E-02	;
//	WV[6].w1[2][8] = 	2.03E+02	;
//	WV[6].w1[2][9] = 	-9.82E-01	;
//	WV[6].w1[2][10] =	1.92E+02	;
//	WV[6].w1[2][11] =	5.48E-01	;
//	WV[6].w1[2][12] =	-2.80E-01	;
//	WV[6].w1[2][13] =	-8.65E+00	;
//	WV[6].w1[2][14] =	3.60E-01	;
//	WV[6].w1[3][0] = 	-6.39E-01	;
//	WV[6].w1[3][1] = 	-2.10E+02	;
//	WV[6].w1[3][2] = 	-1.01E+00	;
//	WV[6].w1[3][3] = 	5.74E-01	;
//	WV[6].w1[3][4] = 	3.72E+01	;
//	WV[6].w1[3][5] = 	-4.36E+00	;
//	WV[6].w1[3][6] = 	5.86E-01	;
//	WV[6].w1[3][7] = 	-4.19E-01	;
//	WV[6].w1[3][8] = 	2.71E+01	;
//	WV[6].w1[3][9] = 	8.58E+01	;
//	WV[6].w1[3][10] =	2.34E+01	;
//	WV[6].w1[3][11] =	3.36E+00	;
//	WV[6].w1[3][12] =	1.44E+00	;
//	WV[6].w1[3][13] =	-3.52E+01	;
//	WV[6].w1[3][14] =	3.49E+00	;
//	WV[6].w2[0][0] = 	2.62E-01	;
//	WV[6].w2[1][0] = 	2.61E+01	;
//	WV[6].w2[2][0] = 	-9.65E+00	;
//	WV[6].w2[3][0] = 	-1.37E+01	;
//	WV[6].w2[4][0] = 	-1.35E-03	;
//	WV[6].w2[5][0] = 	7.35E-02	;
//	WV[6].w2[6][0] = 	1.11E+01	;
//	WV[6].w2[7][0] = 	-1.06E+00	;
//	WV[6].w2[8][0] = 	-1.12E+00	;
//	WV[6].w2[9][0] = 	2.58E+01	;
//	WV[6].w2[10][0] =	1.12E+00	;
//	WV[6].w2[11][0] =	8.60E+00	;
//	WV[6].w2[12][0] =	-1.09E+01	;
//	WV[6].w2[13][0] =	-4.49E-03	;
//	WV[6].w2[14][0] =	-1.46E+01	;
//	WV[6].w2[15][0] =	9.54E+00	;
//	nu_star_MAX[6] = 	1.00E+05	;
//	nu_star_MIN[6] = 	1.00E-12	;
//	D_star_MAX[6]  = 	1.00E+10	;
//	D_star_MIN[6]  = 	1.00E-10	;
//	GGGMAX[6]      = 	3.16E-01	;
WV[6].omega[0] = 	-2.90E+00	;
WV[6].omega[1] = 	-4.34E+00	;
WV[6].omega[2] = 	1.45E+00	;
WV[6].omega[3] = 	3.91E+00	;
WV[6].omega[4] = 	9.44E+00	;
WV[6].omega[5] = 	3.92E+00	;
WV[6].omega[6] = 	4.24E+00	;
WV[6].omega[7] = 	1.33E+01	;
WV[6].omega[8] = 	-2.10E+00	;
WV[6].omega[9] = 	1.92E+00	;
WV[6].omega[10] =	-3.57E+00	;
WV[6].omega[11] =	9.45E+00	;
WV[6].omega[12] =	3.07E+00	;
WV[6].omega[13] =	3.83E+00	;
WV[6].omega[14] =	1.93E+00	;
WV[6].omega[15] =	1.27E+00	;
WV[6].omega[16] =	-7.31E-02	;
WV[6].omega[17] =	1.15E+00	;
WV[6].omega[18] =	-8.52E-01	;
WV[6].omega[19] =	6.30E-01	;
WV[6].omega[20] =	-5.72E-01	;
WV[6].omega[21] =	-5.90E-01	;
WV[6].omega[22] =	-7.69E+00	;
WV[6].omega[23] =	2.22E+00	;
WV[6].omega[24] =	1.19E+00	;
WV[6].omega[25] =	3.73E+00	;
WV[6].omega[26] =	5.78E-01	;
WV[6].omega[27] =	-1.02E+00	;
WV[6].omega[28] =	-4.52E+00	;
WV[6].omega[29] =	-2.76E+00	;
WV[6].omega[30] =	4.05E-01	;
WV[6].omega[31] =	-2.28E-03	;
WV[6].omega[32] =	-9.33E-01	;
WV[6].omega[33] =	-1.31E+00	;
WV[6].omega[34] =	8.51E-01	;
WV[6].omega[35] =	-5.95E-01	;
WV[6].omega[36] =	-8.22E-01	;
WV[6].omega[37] =	-2.00E+00	;
WV[6].omega[38] =	6.86E-01	;
WV[6].omega[39] =	-5.46E-01	;
WV[6].omega[40] =	-1.81E-01	;
WV[6].omega[41] =	8.14E-01	;
WV[6].omega[42] =	-1.57E+00	;
WV[6].omega[43] =	2.40E-01	;
WV[6].omega[44] =	-7.06E-01	;
WV[6].omega[45] =	-4.30E-01	;
WV[6].omega[46] =	1.49E+00	;
WV[6].omega[47] =	2.37E+00	;
WV[6].omega[48] =	1.03E+00	;
WV[6].omega[49] =	4.47E+00	;
WV[6].omega[50] =	8.51E-01	;
WV[6].omega[51] =	9.65E-01	;
WV[6].omega[52] =	9.98E-01	;
WV[6].omega[53] =	-1.83E-01	;
WV[6].omega[54] =	2.08E+00	;
WV[6].omega[55] =	2.27E-01	;
WV[6].omega[56] =	4.43E+00	;
WV[6].omega[57] =	9.69E-01	;
WV[6].omega[58] =	-4.28E-01	;
WV[6].omega[59] =	3.80E-02	;
WV[6].omega[60] =	4.36E+00	;
WV[6].omega[61] =	-3.55E-01	;
WV[6].omega[62] =	1.08E+00	;
WV[6].omega[63] =	2.37E+00	;
WV[6].omega[64] =	-5.00E+00	;
WV[6].omega[65] =	5.90E+00	;
WV[6].omega[66] =	-5.04E+00	;
WV[6].omega[67] =	5.07E-02	;
WV[6].omega[68] =	-3.74E+00	;
WV[6].omega[69] =	-1.49E+00	;
WV[6].omega[70] =	-1.42E+00	;
WV[6].omega[71] =	5.21E+00	;
WV[6].omega[72] =	-1.23E+00	;
WV[6].omega[73] =	-8.74E-01	;
WV[6].omega[74] =	-1.93E+00	;
WV[6].omega[75] =	2.40E-01	;
WV[6].w1[0][0] = 	-2.90E+00	;
WV[6].w1[0][1] = 	-4.34E+00	;
WV[6].w1[0][2] = 	1.45E+00	;
WV[6].w1[0][3] = 	3.91E+00	;
WV[6].w1[0][4] = 	9.44E+00	;
WV[6].w1[0][5] = 	3.92E+00	;
WV[6].w1[0][6] = 	4.24E+00	;
WV[6].w1[0][7] = 	1.33E+01	;
WV[6].w1[0][8] = 	-2.10E+00	;
WV[6].w1[0][9] = 	1.92E+00	;
WV[6].w1[0][10] =	-3.57E+00	;
WV[6].w1[0][11] =	9.45E+00	;
WV[6].w1[0][12] =	3.07E+00	;
WV[6].w1[0][13] =	3.83E+00	;
WV[6].w1[0][14] =	1.93E+00	;
WV[6].w1[1][0] = 	1.27E+00	;
WV[6].w1[1][1] = 	-7.31E-02	;
WV[6].w1[1][2] = 	1.15E+00	;
WV[6].w1[1][3] = 	-8.52E-01	;
WV[6].w1[1][4] = 	6.30E-01	;
WV[6].w1[1][5] = 	-5.72E-01	;
WV[6].w1[1][6] = 	-5.90E-01	;
WV[6].w1[1][7] = 	-7.69E+00	;
WV[6].w1[1][8] = 	2.22E+00	;
WV[6].w1[1][9] = 	1.19E+00	;
WV[6].w1[1][10] =	3.73E+00	;
WV[6].w1[1][11] =	5.78E-01	;
WV[6].w1[1][12] =	-1.02E+00	;
WV[6].w1[1][13] =	-4.52E+00	;
WV[6].w1[1][14] =	-2.76E+00	;
WV[6].w1[2][0] = 	4.05E-01	;
WV[6].w1[2][1] = 	-2.28E-03	;
WV[6].w1[2][2] = 	-9.33E-01	;
WV[6].w1[2][3] = 	-1.31E+00	;
WV[6].w1[2][4] = 	8.51E-01	;
WV[6].w1[2][5] = 	-5.95E-01	;
WV[6].w1[2][6] = 	-8.22E-01	;
WV[6].w1[2][7] = 	-2.00E+00	;
WV[6].w1[2][8] = 	6.86E-01	;
WV[6].w1[2][9] = 	-5.46E-01	;
WV[6].w1[2][10] =	-1.81E-01	;
WV[6].w1[2][11] =	8.14E-01	;
WV[6].w1[2][12] =	-1.57E+00	;
WV[6].w1[2][13] =	2.40E-01	;
WV[6].w1[2][14] =	-7.06E-01	;
WV[6].w1[3][0] = 	-4.30E-01	;
WV[6].w1[3][1] = 	1.49E+00	;
WV[6].w1[3][2] = 	2.37E+00	;
WV[6].w1[3][3] = 	1.03E+00	;
WV[6].w1[3][4] = 	4.47E+00	;
WV[6].w1[3][5] = 	8.51E-01	;
WV[6].w1[3][6] = 	9.65E-01	;
WV[6].w1[3][7] = 	9.98E-01	;
WV[6].w1[3][8] = 	-1.83E-01	;
WV[6].w1[3][9] = 	2.08E+00	;
WV[6].w1[3][10] =	2.27E-01	;
WV[6].w1[3][11] =	4.43E+00	;
WV[6].w1[3][12] =	9.69E-01	;
WV[6].w1[3][13] =	-4.28E-01	;
WV[6].w1[3][14] =	3.80E-02	;
WV[6].w2[0][0] = 	4.36E+00	;
WV[6].w2[1][0] = 	-3.55E-01	;
WV[6].w2[2][0] = 	1.08E+00	;
WV[6].w2[3][0] = 	2.37E+00	;
WV[6].w2[4][0] = 	-5.00E+00	;
WV[6].w2[5][0] = 	5.90E+00	;
WV[6].w2[6][0] = 	-5.04E+00	;
WV[6].w2[7][0] = 	5.07E-02	;
WV[6].w2[8][0] = 	-3.74E+00	;
WV[6].w2[9][0] = 	-1.49E+00	;
WV[6].w2[10][0] =	-1.42E+00	;
WV[6].w2[11][0] =	5.21E+00	;
WV[6].w2[12][0] =	-1.23E+00	;
WV[6].w2[13][0] =	-8.74E-01	;
WV[6].w2[14][0] =	-1.93E+00	;
WV[6].w2[15][0] =	2.40E-01	;
nu_star_MAX[6] = 	1.00E+04	;
nu_star_MIN[6] = 	1.00E-07	;
D_star_MAX[6]  = 	1.00E+05	;
D_star_MIN[6]  = 	1.00E-04	;
GGGMAX[6]      = 	3.16E-01	;




//	/*////////////////   RAX 360 beta 1% data ////////////////*/			
//		WV[6].omega[0] = 	-2.53E+00	;
//		WV[6].omega[1] = 	-3.44E+00	;
//		WV[6].omega[2] = 	-6.90E-01	;
//		WV[6].omega[3] = 	1.65E+00	;
//		WV[6].omega[4] = 	-2.46E+00	;
//		WV[6].omega[5] = 	-3.08E+00	;
//		WV[6].omega[6] = 	-4.49E+00	;
//		WV[6].omega[7] = 	7.73E-02	;
//		WV[6].omega[8] = 	-9.34E+00	;
//		WV[6].omega[9] = 	-1.32E+00	;
//		WV[6].omega[10] =	1.70E+00	;
//		WV[6].omega[11] =	2.50E+00	;
//		WV[6].omega[12] =	4.27E-01	;
//		WV[6].omega[13] =	8.84E+00	;
//		WV[6].omega[14] =	-1.59E+00	;
//		WV[6].omega[15] =	9.03E-01	;
//		WV[6].omega[16] =	1.85E+00	;
//		WV[6].omega[17] =	-4.63E-01	;
//		WV[6].omega[18] =	-2.62E-01	;
//		WV[6].omega[19] =	2.51E+00	;
//		WV[6].omega[20] =	1.96E+00	;
//		WV[6].omega[21] =	3.17E+00	;
//		WV[6].omega[22] =	-6.59E-01	;
//		WV[6].omega[23] =	1.33E+01	;
//		WV[6].omega[24] =	-6.39E-02	;
//		WV[6].omega[25] =	-7.41E-02	;
//		WV[6].omega[26] =	-3.05E+00	;
//		WV[6].omega[27] =	4.95E-01	;
//		WV[6].omega[28] =	-1.38E+00	;
//		WV[6].omega[29] =	3.47E-02	;
//		WV[6].omega[30] =	1.08E+00	;
//		WV[6].omega[31] =	3.72E-01	;
//		WV[6].omega[32] =	-1.65E+00	;
//		WV[6].omega[33] =	-5.25E-01	;
//		WV[6].omega[34] =	4.87E-01	;
//		WV[6].omega[35] =	5.39E-01	;
//		WV[6].omega[36] =	1.11E+00	;
//		WV[6].omega[37] =	3.65E-02	;
//		WV[6].omega[38] =	-4.29E+00	;
//		WV[6].omega[39] =	1.53E+00	;
//		WV[6].omega[40] =	-9.27E-01	;
//		WV[6].omega[41] =	-5.63E-01	;
//		WV[6].omega[42] =	-2.14E-02	;
//		WV[6].omega[43] =	2.34E+01	;
//		WV[6].omega[44] =	4.02E+00	;
//		WV[6].omega[45] =	-7.02E-01	;
//		WV[6].omega[46] =	-5.01E-01	;
//		WV[6].omega[47] =	-2.98E-01	;
//		WV[6].omega[48] =	6.84E-03	;
//		WV[6].omega[49] =	-2.75E-01	;
//		WV[6].omega[50] =	-4.20E-01	;
//		WV[6].omega[51] =	-4.90E-01	;
//		WV[6].omega[52] =	-3.26E-01	;
//		WV[6].omega[53] =	3.85E+00	;
//		WV[6].omega[54] =	-2.61E-01	;
//		WV[6].omega[55] =	4.79E-03	;
//		WV[6].omega[56] =	2.47E-01	;
//		WV[6].omega[57] =	-2.01E-01	;
//		WV[6].omega[58] =	-6.07E+00	;
//		WV[6].omega[59] =	-4.74E-01	;
//		WV[6].omega[60] =	-1.22E+00	;
//		WV[6].omega[61] =	-2.40E+00	;
//		WV[6].omega[62] =	1.30E-01	;
//		WV[6].omega[63] =	-3.12E+00	;
//		WV[6].omega[64] =	-5.08E+00	;
//		WV[6].omega[65] =	6.31E+00	;
//		WV[6].omega[66] =	-6.80E-01	;
//		WV[6].omega[67] =	2.35E+00	;
//		WV[6].omega[68] =	-9.97E-03	;
//		WV[6].omega[69] =	1.55E+00	;
//		WV[6].omega[70] =	2.98E+00	;
//		WV[6].omega[71] =	-2.34E+00	;
//		WV[6].omega[72] =	3.13E+00	;
//		WV[6].omega[73] =	4.15E-02	;
//		WV[6].omega[74] =	-3.12E-01	;
//		WV[6].omega[75] =	1.17E+00	;
//		WV[6].w1[0][0] = 	-2.53E+00	;
//		WV[6].w1[0][1] = 	-3.44E+00	;
//		WV[6].w1[0][2] = 	-6.90E-01	;
//		WV[6].w1[0][3] = 	1.65E+00	;
//		WV[6].w1[0][4] = 	-2.46E+00	;
//		WV[6].w1[0][5] = 	-3.08E+00	;
//		WV[6].w1[0][6] = 	-4.49E+00	;
//		WV[6].w1[0][7] = 	7.73E-02	;
//		WV[6].w1[0][8] = 	-9.34E+00	;
//		WV[6].w1[0][9] = 	-1.32E+00	;
//		WV[6].w1[0][10] =	1.70E+00	;
//		WV[6].w1[0][11] =	2.50E+00	;
//		WV[6].w1[0][12] =	4.27E-01	;
//		WV[6].w1[0][13] =	8.84E+00	;
//		WV[6].w1[0][14] =	-1.59E+00	;
//		WV[6].w1[1][0] = 	9.03E-01	;
//		WV[6].w1[1][1] = 	1.85E+00	;
//		WV[6].w1[1][2] = 	-4.63E-01	;
//		WV[6].w1[1][3] = 	-2.62E-01	;
//		WV[6].w1[1][4] = 	2.51E+00	;
//		WV[6].w1[1][5] = 	1.96E+00	;
//		WV[6].w1[1][6] = 	3.17E+00	;
//		WV[6].w1[1][7] = 	-6.59E-01	;
//		WV[6].w1[1][8] = 	1.33E+01	;
//		WV[6].w1[1][9] = 	-6.39E-02	;
//		WV[6].w1[1][10] =	-7.41E-02	;
//		WV[6].w1[1][11] =	-3.05E+00	;
//		WV[6].w1[1][12] =	4.95E-01	;
//		WV[6].w1[1][13] =	-1.38E+00	;
//		WV[6].w1[1][14] =	3.47E-02	;
//		WV[6].w1[2][0] = 	1.08E+00	;
//		WV[6].w1[2][1] = 	3.72E-01	;
//		WV[6].w1[2][2] = 	-1.65E+00	;
//		WV[6].w1[2][3] = 	-5.25E-01	;
//		WV[6].w1[2][4] = 	4.87E-01	;
//		WV[6].w1[2][5] = 	5.39E-01	;
//		WV[6].w1[2][6] = 	1.11E+00	;
//		WV[6].w1[2][7] = 	3.65E-02	;
//		WV[6].w1[2][8] = 	-4.29E+00	;
//		WV[6].w1[2][9] = 	1.53E+00	;
//		WV[6].w1[2][10] =	-9.27E-01	;
//		WV[6].w1[2][11] =	-5.63E-01	;
//		WV[6].w1[2][12] =	-2.14E-02	;
//		WV[6].w1[2][13] =	2.34E+01	;
//		WV[6].w1[2][14] =	4.02E+00	;
//		WV[6].w1[3][0] = 	-7.02E-01	;
//		WV[6].w1[3][1] = 	-5.01E-01	;
//		WV[6].w1[3][2] = 	-2.98E-01	;
//		WV[6].w1[3][3] = 	6.84E-03	;
//		WV[6].w1[3][4] = 	-2.75E-01	;
//		WV[6].w1[3][5] = 	-4.20E-01	;
//		WV[6].w1[3][6] = 	-4.90E-01	;
//		WV[6].w1[3][7] = 	-3.26E-01	;
//		WV[6].w1[3][8] = 	3.85E+00	;
//		WV[6].w1[3][9] = 	-2.61E-01	;
//		WV[6].w1[3][10] =	4.79E-03	;
//		WV[6].w1[3][11] =	2.47E-01	;
//		WV[6].w1[3][12] =	-2.01E-01	;
//		WV[6].w1[3][13] =	-6.07E+00	;
//		WV[6].w1[3][14] =	-4.74E-01	;
//		WV[6].w2[0][0] = 	-1.22E+00	;
//		WV[6].w2[1][0] = 	-2.40E+00	;
//		WV[6].w2[2][0] = 	1.30E-01	;
//		WV[6].w2[3][0] = 	-3.12E+00	;
//		WV[6].w2[4][0] = 	-5.08E+00	;
//		WV[6].w2[5][0] = 	6.31E+00	;
//		WV[6].w2[6][0] = 	-6.80E-01	;
//		WV[6].w2[7][0] = 	2.35E+00	;
//		WV[6].w2[8][0] = 	-9.97E-03	;
//		WV[6].w2[9][0] = 	1.55E+00	;
//		WV[6].w2[10][0] =	2.98E+00	;
//		WV[6].w2[11][0] =	-2.34E+00	;
//		WV[6].w2[12][0] =	3.13E+00	;
//		WV[6].w2[13][0] =	4.15E-02	;
//		WV[6].w2[14][0] =	-3.12E-01	;
//		WV[6].w2[15][0] =	1.17E+00	;
//		nu_star_MAX[6] = 	1.00E+03	;
//		nu_star_MIN[6] = 	1.00E-06	;
//		D_star_MAX[6]  = 	1.00E+05	;
//		D_star_MIN[6]  = 	1.00E-04	;
//		GGGMAX[6]      = 	3.16E-01	;


/*////////////////   RAX 360 beta 2% data ////////////////*/
WV[7].omega[0] = 	6.99748E-01	;
WV[7].omega[1] = 	1.76662E+00	;
WV[7].omega[2] = 	-1.20348E+00	;
WV[7].omega[3] = 	-3.64275E+00	;
WV[7].omega[4] = 	-1.36135E+00	;
WV[7].omega[5] = 	1.19071E+00	;
WV[7].omega[6] = 	1.07006E+00	;
WV[7].omega[7] = 	-4.11813E+00	;
WV[7].omega[8] = 	-4.26836E+00	;
WV[7].omega[9] = 	9.68198E-01	;
WV[7].omega[10] =	-4.75705E+00	;
WV[7].omega[11] =	-3.84815E+00	;
WV[7].omega[12] =	-4.03984E+00	;
WV[7].omega[13] =	-1.35555E+00	;
WV[7].omega[14] =	-1.17866E+00	;
WV[7].omega[15] =	-3.60760E-01	;
WV[7].omega[16] =	-3.26830E+00	;
WV[7].omega[17] =	-3.55403E+00	;
WV[7].omega[18] =	-5.36085E-02	;
WV[7].omega[19] =	9.04523E-01	;
WV[7].omega[20] =	6.79239E+00	;
WV[7].omega[21] =	-1.12440E+00	;
WV[7].omega[22] =	5.67212E+00	;
WV[7].omega[23] =	3.68918E+00	;
WV[7].omega[24] =	-1.34643E+00	;
WV[7].omega[25] =	3.63351E+00	;
WV[7].omega[26] =	-2.39998E-02	;
WV[7].omega[27] =	4.80423E+00	;
WV[7].omega[28] =	1.45499E-01	;
WV[7].omega[29] =	9.61831E-01	;
WV[7].omega[30] =	4.11923E+00	;
WV[7].omega[31] =	-1.86623E-01	;
WV[7].omega[32] =	8.53342E-02	;
WV[7].omega[33] =	-1.48107E-01	;
WV[7].omega[34] =	7.68835E-02	;
WV[7].omega[35] =	-5.14607E-02	;
WV[7].omega[36] =	6.83744E-01	;
WV[7].omega[37] =	5.19618E-01	;
WV[7].omega[38] =	4.83733E-01	;
WV[7].omega[39] =	1.25502E+00	;
WV[7].omega[40] =	5.32511E-01	;
WV[7].omega[41] =	5.04814E-02	;
WV[7].omega[42] =	4.94745E-01	;
WV[7].omega[43] =	-2.44271E+00	;
WV[7].omega[44] =	2.97055E-01	;
WV[7].omega[45] =	-2.16950E+00	;
WV[7].omega[46] =	1.17125E+00	;
WV[7].omega[47] =	-4.74990E+00	;
WV[7].omega[48] =	4.61058E+00	;
WV[7].omega[49] =	1.67181E-02	;
WV[7].omega[50] =	7.66747E+00	;
WV[7].omega[51] =	4.96075E-01	;
WV[7].omega[52] =	6.39022E-01	;
WV[7].omega[53] =	2.05167E-01	;
WV[7].omega[54] =	6.85387E-01	;
WV[7].omega[55] =	5.30293E-02	;
WV[7].omega[56] =	1.52061E+00	;
WV[7].omega[57] =	4.90900E-01	;
WV[7].omega[58] =	1.28399E+00	;
WV[7].omega[59] =	-5.67068E-01	;
WV[7].omega[60] =	4.92473E-02	;
WV[7].omega[61] =	-4.93436E-01	;
WV[7].omega[62] =	3.65039E+00	;
WV[7].omega[63] =	-2.50879E+00	;
WV[7].omega[64] =	1.37871E+00	;
WV[7].omega[65] =	4.26607E+00	;
WV[7].omega[66] =	1.33871E+00	;
WV[7].omega[67] =	7.25591E-01	;
WV[7].omega[68] =	1.96417E+00	;
WV[7].omega[69] =	-5.32677E-01	;
WV[7].omega[70] =	-1.06546E+00	;
WV[7].omega[71] =	-2.86851E-01	;
WV[7].omega[72] =	-1.56253E+00	;
WV[7].omega[73] =	4.40680E-02	;
WV[7].omega[74] =	-1.87272E+00	;
WV[7].omega[75] =	1.18003E+00	;
WV[7].w1[0][0] = 	6.99748E-01	;
WV[7].w1[0][1] = 	1.76662E+00	;
WV[7].w1[0][2] = 	-1.20348E+00	;
WV[7].w1[0][3] = 	-3.64275E+00	;
WV[7].w1[0][4] = 	-1.36135E+00	;
WV[7].w1[0][5] = 	1.19071E+00	;
WV[7].w1[0][6] = 	1.07006E+00	;
WV[7].w1[0][7] = 	-4.11813E+00	;
WV[7].w1[0][8] = 	-4.26836E+00	;
WV[7].w1[0][9] = 	9.68198E-01	;
WV[7].w1[0][10] =	-4.75705E+00	;
WV[7].w1[0][11] =	-3.84815E+00	;
WV[7].w1[0][12] =	-4.03984E+00	;
WV[7].w1[0][13] =	-1.35555E+00	;
WV[7].w1[0][14] =	-1.17866E+00	;
WV[7].w1[1][0] = 	-3.60760E-01	;
WV[7].w1[1][1] = 	-3.26830E+00	;
WV[7].w1[1][2] = 	-3.55403E+00	;
WV[7].w1[1][3] = 	-5.36085E-02	;
WV[7].w1[1][4] = 	9.04523E-01	;
WV[7].w1[1][5] = 	6.79239E+00	;
WV[7].w1[1][6] = 	-1.12440E+00	;
WV[7].w1[1][7] = 	5.67212E+00	;
WV[7].w1[1][8] = 	3.68918E+00	;
WV[7].w1[1][9] = 	-1.34643E+00	;
WV[7].w1[1][10] =	3.63351E+00	;
WV[7].w1[1][11] =	-2.39998E-02	;
WV[7].w1[1][12] =	4.80423E+00	;
WV[7].w1[1][13] =	1.45499E-01	;
WV[7].w1[1][14] =	9.61831E-01	;
WV[7].w1[2][0] = 	4.11923E+00	;
WV[7].w1[2][1] = 	-1.86623E-01	;
WV[7].w1[2][2] = 	8.53342E-02	;
WV[7].w1[2][3] = 	-1.48107E-01	;
WV[7].w1[2][4] = 	7.68835E-02	;
WV[7].w1[2][5] = 	-5.14607E-02	;
WV[7].w1[2][6] = 	6.83744E-01	;
WV[7].w1[2][7] = 	5.19618E-01	;
WV[7].w1[2][8] = 	4.83733E-01	;
WV[7].w1[2][9] = 	1.25502E+00	;
WV[7].w1[2][10] =	5.32511E-01	;
WV[7].w1[2][11] =	5.04814E-02	;
WV[7].w1[2][12] =	4.94745E-01	;
WV[7].w1[2][13] =	-2.44271E+00	;
WV[7].w1[2][14] =	2.97055E-01	;
WV[7].w1[3][0] = 	-2.16950E+00	;
WV[7].w1[3][1] = 	1.17125E+00	;
WV[7].w1[3][2] = 	-4.74990E+00	;
WV[7].w1[3][3] = 	4.61058E+00	;
WV[7].w1[3][4] = 	1.67181E-02	;
WV[7].w1[3][5] = 	7.66747E+00	;
WV[7].w1[3][6] = 	4.96075E-01	;
WV[7].w1[3][7] = 	6.39022E-01	;
WV[7].w1[3][8] = 	2.05167E-01	;
WV[7].w1[3][9] = 	6.85387E-01	;
WV[7].w1[3][10] =	5.30293E-02	;
WV[7].w1[3][11] =	1.52061E+00	;
WV[7].w1[3][12] =	4.90900E-01	;
WV[7].w1[3][13] =	1.28399E+00	;
WV[7].w1[3][14] =	-5.67068E-01	;
WV[7].w2[0][0] = 	4.92473E-02	;
WV[7].w2[1][0] = 	-4.93436E-01	;
WV[7].w2[2][0] = 	3.65039E+00	;
WV[7].w2[3][0] = 	-2.50879E+00	;
WV[7].w2[4][0] = 	1.37871E+00	;
WV[7].w2[5][0] = 	4.26607E+00	;
WV[7].w2[6][0] = 	1.33871E+00	;
WV[7].w2[7][0] = 	7.25591E-01	;
WV[7].w2[8][0] = 	1.96417E+00	;
WV[7].w2[9][0] = 	-5.32677E-01	;
WV[7].w2[10][0] =	-1.06546E+00	;
WV[7].w2[11][0] =	-2.86851E-01	;
WV[7].w2[12][0] =	-1.56253E+00	;
WV[7].w2[13][0] =	4.40680E-02	;
WV[7].w2[14][0] =	-1.87272E+00	;
WV[7].w2[15][0] =	1.18003E+00	;
nu_star_MAX[7] = 	1.000000E+03	;
nu_star_MIN[7] = 	1.000000E-06	;
D_star_MAX[7]  = 	1.000000E+05	;
D_star_MIN[7]  = 	1.000000E-04	;
GGGMAX[7]      = 	3.160000E-01	;

/*////////////////   RAX 360 beta 3% data ////////////////*/
WV[8].omega[0] = 	-2.038201E+00	;
WV[8].omega[1] = 	3.791647E+00	;
WV[8].omega[2] = 	-1.358189E+00	;
WV[8].omega[3] = 	-4.463407E+00	;
WV[8].omega[4] = 	3.694268E+00	;
WV[8].omega[5] = 	2.289873E+00	;
WV[8].omega[6] = 	8.792521E+00	;
WV[8].omega[7] = 	-5.478823E+00	;
WV[8].omega[8] = 	-3.236900E+00	;
WV[8].omega[9] = 	8.542165E+00	;
WV[8].omega[10] =	3.656371E+00	;
WV[8].omega[11] =	1.973440E+00	;
WV[8].omega[12] =	1.167526E+00	;
WV[8].omega[13] =	2.223781E+00	;
WV[8].omega[14] =	-1.264316E+00	;
WV[8].omega[15] =	-2.852807E+00	;
WV[8].omega[16] =	2.317845E+01	;
WV[8].omega[17] =	-5.453711E-03	;
WV[8].omega[18] =	4.871454E+00	;
WV[8].omega[19] =	-2.964366E+00	;
WV[8].omega[20] =	1.482181E+00	;
WV[8].omega[21] =	1.674791E-02	;
WV[8].omega[22] =	-9.124849E-02	;
WV[8].omega[23] =	-6.007064E-01	;
WV[8].omega[24] =	3.188256E-02	;
WV[8].omega[25] =	-3.002603E+00	;
WV[8].omega[26] =	-3.781530E-02	;
WV[8].omega[27] =	-1.678122E-01	;
WV[8].omega[28] =	6.785658E+00	;
WV[8].omega[29] =	1.061095E+00	;
WV[8].omega[30] =	6.062659E-01	;
WV[8].omega[31] =	-3.822463E-01	;
WV[8].omega[32] =	-7.349980E-02	;
WV[8].omega[33] =	-1.004780E+00	;
WV[8].omega[34] =	-3.531589E-01	;
WV[8].omega[35] =	-5.875387E-01	;
WV[8].omega[36] =	5.535616E+00	;
WV[8].omega[37] =	1.131342E+00	;
WV[8].omega[38] =	6.409435E-01	;
WV[8].omega[39] =	5.164052E+00	;
WV[8].omega[40] =	-3.479815E-01	;
WV[8].omega[41] =	-1.799154E-02	;
WV[8].omega[42] =	-2.947161E-01	;
WV[8].omega[43] =	-1.805639E-01	;
WV[8].omega[44] =	3.690674E-01	;
WV[8].omega[45] =	-3.179748E+00	;
WV[8].omega[46] =	2.009880E+01	;
WV[8].omega[47] =	4.660970E-01	;
WV[8].omega[48] =	-6.920744E-01	;
WV[8].omega[49] =	-7.007160E-01	;
WV[8].omega[50] =	2.652907E+00	;
WV[8].omega[51] =	-6.635768E+00	;
WV[8].omega[52] =	4.788580E-01	;
WV[8].omega[53] =	-1.218771E+00	;
WV[8].omega[54] =	-6.312176E+00	;
WV[8].omega[55] =	-6.897645E-01	;
WV[8].omega[56] =	-3.353893E-01	;
WV[8].omega[57] =	1.311956E+00	;
WV[8].omega[58] =	7.891964E+00	;
WV[8].omega[59] =	-2.044178E+00	;
WV[8].omega[60] =	-7.144257E-01	;
WV[8].omega[61] =	3.386478E-01	;
WV[8].omega[62] =	-1.629490E+00	;
WV[8].omega[63] =	-4.745556E-02	;
WV[8].omega[64] =	-3.962794E+00	;
WV[8].omega[65] =	-2.310204E+00	;
WV[8].omega[66] =	7.087521E-01	;
WV[8].omega[67] =	8.628787E-02	;
WV[8].omega[68] =	4.134925E-01	;
WV[8].omega[69] =	-7.472096E-01	;
WV[8].omega[70] =	3.950671E+00	;
WV[8].omega[71] =	-1.095076E+00	;
WV[8].omega[72] =	6.157182E+00	;
WV[8].omega[73] =	-3.704954E+00	;
WV[8].omega[74] =	6.005810E+00	;
WV[8].omega[75] =	5.960201E+00	;
WV[8].w1[0][0] = 	-2.038201E+00	;
WV[8].w1[0][1] = 	3.791647E+00	;
WV[8].w1[0][2] = 	-1.358189E+00	;
WV[8].w1[0][3] = 	-4.463407E+00	;
WV[8].w1[0][4] = 	3.694268E+00	;
WV[8].w1[0][5] = 	2.289873E+00	;
WV[8].w1[0][6] = 	8.792521E+00	;
WV[8].w1[0][7] = 	-5.478823E+00	;
WV[8].w1[0][8] = 	-3.236900E+00	;
WV[8].w1[0][9] = 	8.542165E+00	;
WV[8].w1[0][10] =	3.656371E+00	;
WV[8].w1[0][11] =	1.973440E+00	;
WV[8].w1[0][12] =	1.167526E+00	;
WV[8].w1[0][13] =	2.223781E+00	;
WV[8].w1[0][14] =	-1.264316E+00	;
WV[8].w1[1][0] = 	-2.852807E+00	;
WV[8].w1[1][1] = 	2.317845E+01	;
WV[8].w1[1][2] = 	-5.453711E-03	;
WV[8].w1[1][3] = 	4.871454E+00	;
WV[8].w1[1][4] = 	-2.964366E+00	;
WV[8].w1[1][5] = 	1.482181E+00	;
WV[8].w1[1][6] = 	1.674791E-02	;
WV[8].w1[1][7] = 	-9.124849E-02	;
WV[8].w1[1][8] = 	-6.007064E-01	;
WV[8].w1[1][9] = 	3.188256E-02	;
WV[8].w1[1][10] =	-3.002603E+00	;
WV[8].w1[1][11] =	-3.781530E-02	;
WV[8].w1[1][12] =	-1.678122E-01	;
WV[8].w1[1][13] =	6.785658E+00	;
WV[8].w1[1][14] =	1.061095E+00	;
WV[8].w1[2][0] = 	6.062659E-01	;
WV[8].w1[2][1] = 	-3.822463E-01	;
WV[8].w1[2][2] = 	-7.349980E-02	;
WV[8].w1[2][3] = 	-1.004780E+00	;
WV[8].w1[2][4] = 	-3.531589E-01	;
WV[8].w1[2][5] = 	-5.875387E-01	;
WV[8].w1[2][6] = 	5.535616E+00	;
WV[8].w1[2][7] = 	1.131342E+00	;
WV[8].w1[2][8] = 	6.409435E-01	;
WV[8].w1[2][9] = 	5.164052E+00	;
WV[8].w1[2][10] =	-3.479815E-01	;
WV[8].w1[2][11] =	-1.799154E-02	;
WV[8].w1[2][12] =	-2.947161E-01	;
WV[8].w1[2][13] =	-1.805639E-01	;
WV[8].w1[2][14] =	3.690674E-01	;
WV[8].w1[3][0] = 	-3.179748E+00	;
WV[8].w1[3][1] = 	2.009880E+01	;
WV[8].w1[3][2] = 	4.660970E-01	;
WV[8].w1[3][3] = 	-6.920744E-01	;
WV[8].w1[3][4] = 	-7.007160E-01	;
WV[8].w1[3][5] = 	2.652907E+00	;
WV[8].w1[3][6] = 	-6.635768E+00	;
WV[8].w1[3][7] = 	4.788580E-01	;
WV[8].w1[3][8] = 	-1.218771E+00	;
WV[8].w1[3][9] = 	-6.312176E+00	;
WV[8].w1[3][10] =	-6.897645E-01	;
WV[8].w1[3][11] =	-3.353893E-01	;
WV[8].w1[3][12] =	1.311956E+00	;
WV[8].w1[3][13] =	7.891964E+00	;
WV[8].w1[3][14] =	-2.044178E+00	;
WV[8].w2[0][0] = 	-7.144257E-01	;
WV[8].w2[1][0] = 	3.386478E-01	;
WV[8].w2[2][0] = 	-1.629490E+00	;
WV[8].w2[3][0] = 	-4.745556E-02	;
WV[8].w2[4][0] = 	-3.962794E+00	;
WV[8].w2[5][0] = 	-2.310204E+00	;
WV[8].w2[6][0] = 	7.087521E-01	;
WV[8].w2[7][0] = 	8.628787E-02	;
WV[8].w2[8][0] = 	4.134925E-01	;
WV[8].w2[9][0] = 	-7.472096E-01	;
WV[8].w2[10][0] =	3.950671E+00	;
WV[8].w2[11][0] =	-1.095076E+00	;
WV[8].w2[12][0] =	6.157182E+00	;
WV[8].w2[13][0] =	-3.704954E+00	;
WV[8].w2[14][0] =	6.005810E+00	;
WV[8].w2[15][0] =	5.960201E+00	;
nu_star_MAX[8] = 	1.00000E+03	;
nu_star_MIN[8] = 	1.00000E-06	;
D_star_MAX[8]  = 	1.00000E+05	;
D_star_MIN[8]  = 	1.00000E-04	;
GGGMAX[8]      = 	3.160000E-01	;

/*////////////////   RAX 360 beta 0.5% data ////////////////*/		
WV[9].omega[0] = 	2.937418E+00	;
WV[9].omega[1] = 	3.285380E+00	;
WV[9].omega[2] = 	3.657868E+00	;
WV[9].omega[3] = 	-4.361186E+00	;
WV[9].omega[4] = 	6.221446E+00	;
WV[9].omega[5] = 	-1.299542E+00	;
WV[9].omega[6] = 	1.658842E+01	;
WV[9].omega[7] = 	-2.514283E+00	;
WV[9].omega[8] = 	3.862041E+00	;
WV[9].omega[9] = 	4.321528E+00	;
WV[9].omega[10] =	7.001490E+00	;
WV[9].omega[11] =	2.973911E+00	;
WV[9].omega[12] =	-2.861459E+00	;
WV[9].omega[13] =	-2.972982E+00	;
WV[9].omega[14] =	-3.162675E+00	;
WV[9].omega[15] =	-4.982707E+00	;
WV[9].omega[16] =	-7.946919E-01	;
WV[9].omega[17] =	-1.545923E+00	;
WV[9].omega[18] =	8.959147E-01	;
WV[9].omega[19] =	1.301760E-01	;
WV[9].omega[20] =	1.343960E+00	;
WV[9].omega[21] =	-2.517946E-02	;
WV[9].omega[22] =	1.829383E+00	;
WV[9].omega[23] =	2.512452E+00	;
WV[9].omega[24] =	-7.787690E-01	;
WV[9].omega[25] =	2.161035E-01	;
WV[9].omega[26] =	-4.214760E+00	;
WV[9].omega[27] =	-3.072870E-01	;
WV[9].omega[28] =	5.784523E+00	;
WV[9].omega[29] =	4.940800E+00	;
WV[9].omega[30] =	-4.373775E-01	;
WV[9].omega[31] =	9.396767E-01	;
WV[9].omega[32] =	1.688554E+00	;
WV[9].omega[33] =	1.535407E+00	;
WV[9].omega[34] =	-1.963413E-01	;
WV[9].omega[35] =	2.288664E-01	;
WV[9].omega[36] =	5.933072E-01	;
WV[9].omega[37] =	3.295896E-01	;
WV[9].omega[38] =	-2.914352E-01	;
WV[9].omega[39] =	-1.435746E+00	;
WV[9].omega[40] =	-1.168111E+00	;
WV[9].omega[41] =	-4.361871E-01	;
WV[9].omega[42] =	9.810159E-06	;
WV[9].omega[43] =	4.673638E-01	;
WV[9].omega[44] =	4.808780E-01	;
WV[9].omega[45] =	-1.771276E+00	;
WV[9].omega[46] =	2.683552E+00	;
WV[9].omega[47] =	2.696377E+00	;
WV[9].omega[48] =	-6.178078E-01	;
WV[9].omega[49] =	-3.098958E+00	;
WV[9].omega[50] =	4.774535E-01	;
WV[9].omega[51] =	-4.027113E+00	;
WV[9].omega[52] =	1.491915E+00	;
WV[9].omega[53] =	3.652687E+00	;
WV[9].omega[54] =	6.204337E-01	;
WV[9].omega[55] =	9.238757E-01	;
WV[9].omega[56] =	-1.281856E+00	;
WV[9].omega[57] =	1.565927E+00	;
WV[9].omega[58] =	1.932128E+00	;
WV[9].omega[59] =	1.426613E+00	;
WV[9].omega[60] =	2.202763E+00	;
WV[9].omega[61] =	2.652053E+00	;
WV[9].omega[62] =	-9.420637E-01	;
WV[9].omega[63] =	-1.165089E+00	;
WV[9].omega[64] =	-2.880455E-01	;
WV[9].omega[65] =	-1.622309E+00	;
WV[9].omega[66] =	2.392592E-02	;
WV[9].omega[67] =	6.008991E-01	;
WV[9].omega[68] =	-1.610062E-01	;
WV[9].omega[69] =	-1.383241E+00	;
WV[9].omega[70] =	9.768143E-02	;
WV[9].omega[71] =	-4.580348E+00	;
WV[9].omega[72] =	-1.200141E+00	;
WV[9].omega[73] =	1.431643E+00	;
WV[9].omega[74] =	-3.159820E+00	;
WV[9].omega[75] =	-9.794720E-01	;
WV[9].w1[0][0] = 	2.937418E+00	;
WV[9].w1[0][1] = 	3.285380E+00	;
WV[9].w1[0][2] = 	3.657868E+00	;
WV[9].w1[0][3] = 	-4.361186E+00	;
WV[9].w1[0][4] = 	6.221446E+00	;
WV[9].w1[0][5] = 	-1.299542E+00	;
WV[9].w1[0][6] = 	1.658842E+01	;
WV[9].w1[0][7] = 	-2.514283E+00	;
WV[9].w1[0][8] = 	3.862041E+00	;
WV[9].w1[0][9] = 	4.321528E+00	;
WV[9].w1[0][10] =	7.001490E+00	;
WV[9].w1[0][11] =	2.973911E+00	;
WV[9].w1[0][12] =	-2.861459E+00	;
WV[9].w1[0][13] =	-2.972982E+00	;
WV[9].w1[0][14] =	-3.162675E+00	;
WV[9].w1[1][0] = 	-4.982707E+00	;
WV[9].w1[1][1] = 	-7.946919E-01	;
WV[9].w1[1][2] = 	-1.545923E+00	;
WV[9].w1[1][3] = 	8.959147E-01	;
WV[9].w1[1][4] = 	1.301760E-01	;
WV[9].w1[1][5] = 	1.343960E+00	;
WV[9].w1[1][6] = 	-2.517946E-02	;
WV[9].w1[1][7] = 	1.829383E+00	;
WV[9].w1[1][8] = 	2.512452E+00	;
WV[9].w1[1][9] = 	-7.787690E-01	;
WV[9].w1[1][10] =	2.161035E-01	;
WV[9].w1[1][11] =	-4.214760E+00	;
WV[9].w1[1][12] =	-3.072870E-01	;
WV[9].w1[1][13] =	5.784523E+00	;
WV[9].w1[1][14] =	4.940800E+00	;
WV[9].w1[2][0] = 	-4.373775E-01	;
WV[9].w1[2][1] = 	9.396767E-01	;
WV[9].w1[2][2] = 	1.688554E+00	;
WV[9].w1[2][3] = 	1.535407E+00	;
WV[9].w1[2][4] = 	-1.963413E-01	;
WV[9].w1[2][5] = 	2.288664E-01	;
WV[9].w1[2][6] = 	5.933072E-01	;
WV[9].w1[2][7] = 	3.295896E-01	;
WV[9].w1[2][8] = 	-2.914352E-01	;
WV[9].w1[2][9] = 	-1.435746E+00	;
WV[9].w1[2][10] =	-1.168111E+00	;
WV[9].w1[2][11] =	-4.361871E-01	;
WV[9].w1[2][12] =	9.810159E-06	;
WV[9].w1[2][13] =	4.673638E-01	;
WV[9].w1[2][14] =	4.808780E-01	;
WV[9].w1[3][0] = 	-1.771276E+00	;
WV[9].w1[3][1] = 	2.683552E+00	;
WV[9].w1[3][2] = 	2.696377E+00	;
WV[9].w1[3][3] = 	-6.178078E-01	;
WV[9].w1[3][4] = 	-3.098958E+00	;
WV[9].w1[3][5] = 	4.774535E-01	;
WV[9].w1[3][6] = 	-4.027113E+00	;
WV[9].w1[3][7] = 	1.491915E+00	;
WV[9].w1[3][8] = 	3.652687E+00	;
WV[9].w1[3][9] = 	6.204337E-01	;
WV[9].w1[3][10] =	9.238757E-01	;
WV[9].w1[3][11] =	-1.281856E+00	;
WV[9].w1[3][12] =	1.565927E+00	;
WV[9].w1[3][13] =	1.932128E+00	;
WV[9].w1[3][14] =	1.426613E+00	;
WV[9].w2[0][0] = 	2.202763E+00	;
WV[9].w2[1][0] = 	2.652053E+00	;
WV[9].w2[2][0] = 	-9.420637E-01	;
WV[9].w2[3][0] = 	-1.165089E+00	;
WV[9].w2[4][0] = 	-2.880455E-01	;
WV[9].w2[5][0] = 	-1.622309E+00	;
WV[9].w2[6][0] = 	2.392592E-02	;
WV[9].w2[7][0] = 	6.008991E-01	;
WV[9].w2[8][0] = 	-1.610062E-01	;
WV[9].w2[9][0] = 	-1.383241E+00	;
WV[9].w2[10][0] =	9.768143E-02	;
WV[9].w2[11][0] =	-4.580348E+00	;
WV[9].w2[12][0] =	-1.200141E+00	;
WV[9].w2[13][0] =	1.431643E+00	;
WV[9].w2[14][0] =	-3.159820E+00	;
WV[9].w2[15][0] =	-9.794720E-01	;
nu_star_MAX[9] = 	1.00E+04	;
nu_star_MIN[9] = 	1.00E-07	;
D_star_MAX[9]  = 	1.00E+05	;
D_star_MIN[9]  = 	1.00E-04	;
GGGMAX[9]      = 	3.16E-01	;

#ifdef use_NNWwithBETA
/*////////////////   RAX 375 finite beta data ////////////////*/
WV[15].omega[0] =	-1.728066E+01	;
WV[15].omega[1] =	2.210435E+00	;
WV[15].omega[2] =	2.348677E+00	;
WV[15].omega[3] =	2.812618E-01	;
WV[15].omega[4] =	2.417725E+00	;
WV[15].omega[5] =	1.246956E+00	;
WV[15].omega[6] =	-6.094876E-01	;
WV[15].omega[7] =	-5.309743E+00	;
WV[15].omega[8] =	-2.152382E+00	;
WV[15].omega[9] =	1.873751E+00	;
WV[15].omega[10] =	7.651818E-02	;
WV[15].omega[11] =	-1.529941E+01	;
WV[15].omega[12] =	1.839610E+00	;
WV[15].omega[13] =	-2.523251E+00	;
WV[15].omega[14] =	2.428332E+00	;
WV[15].omega[15] =	-2.273623E-02	;
WV[15].omega[16] =	-1.737884E+00	;
WV[15].omega[17] =	7.905896E-01	;
WV[15].omega[18] =	6.848911E-01	;
WV[15].omega[19] =	-3.599874E+00	;
WV[15].omega[20] =	-6.741247E-02	;
WV[15].omega[21] =	9.020591E-01	;
WV[15].omega[22] =	-7.354806E-03	;
WV[15].omega[23] =	9.320515E-02	;
WV[15].omega[24] =	1.300522E+00	;
WV[15].omega[25] =	-7.692179E-01	;
WV[15].omega[26] =	2.020262E+00	;
WV[15].omega[27] =	1.457440E+00	;
WV[15].omega[28] =	3.024908E+00	;
WV[15].omega[29] =	-4.385509E+00	;
WV[15].omega[30] =	5.227908E-01	;
WV[15].omega[31] =	6.498177E-01	;
WV[15].omega[32] =	-3.313994E-01	;
WV[15].omega[33] =	6.635243E-01	;
WV[15].omega[34] =	-3.528461E-01	;
WV[15].omega[35] =	-3.097976E-01	;
WV[15].omega[36] =	6.269263E-01	;
WV[15].omega[37] =	1.979896E-01	;
WV[15].omega[38] =	-5.230736E-01	;
WV[15].omega[39] =	-5.516489E-01	;
WV[15].omega[40] =	-6.330952E-01	;
WV[15].omega[41] =	1.594426E+01	;
WV[15].omega[42] =	-6.062846E-01	;
WV[15].omega[43] =	3.653671E-01	;
WV[15].omega[44] =	-3.582658E-01	;
WV[15].omega[45] =	2.421305E+00	;
WV[15].omega[46] =	7.068348E-01	;
WV[15].omega[47] =	-2.395018E+00	;
WV[15].omega[48] =	2.980795E+00	;
WV[15].omega[49] =	-1.885558E+00	;
WV[15].omega[50] =	-2.182956E+00	;
WV[15].omega[51] =	3.518668E+00	;
WV[15].omega[52] =	1.327941E+00	;
WV[15].omega[53] =	-2.535836E+00	;
WV[15].omega[54] =	-4.720176E+00	;
WV[15].omega[55] =	-3.051918E+00	;
WV[15].omega[56] =	-6.070693E+00	;
WV[15].omega[57] =	-5.186735E+00	;
WV[15].omega[58] =	2.165606E+00	;
WV[15].omega[59] =	-1.691555E+00	;
WV[15].omega[60] =	1.399478E+01	;
WV[15].omega[61] =	2.389718E+00	;
WV[15].omega[62] =	1.780490E+00	;
WV[15].omega[63] =	2.466100E+00	;
WV[15].omega[64] =	-2.008582E+00	;
WV[15].omega[65] =	2.876683E-01	;
WV[15].omega[66] =	3.439111E+00	;
WV[15].omega[67] =	3.260034E+00	;
WV[15].omega[68] =	-3.673472E+00	;
WV[15].omega[69] =	-1.727469E+00	;
WV[15].omega[70] =	-2.751030E+00	;
WV[15].omega[71] =	2.712121E+00	;
WV[15].omega[72] =	-1.946774E+00	;
WV[15].omega[73] =	2.165539E+00	;
WV[15].omega[74] =	-1.973277E+00	;
WV[15].omega[75] =	-1.152547E-01	;
WV[15].omega[76] =	-4.966613E-01	;
WV[15].omega[77] =	-5.353258E+00	;
WV[15].omega[78] =	1.614083E+00	;
WV[15].omega[79] =	3.080289E+00	;
WV[15].omega[80] =	3.921723E+00	;
WV[15].omega[81] =	1.062526E+00	;
WV[15].omega[82] =	-2.159306E-01	;
WV[15].omega[83] =	-8.224340E-01	;
WV[15].omega[84] =	-1.889153E+00	;
WV[15].omega[85] =	2.619577E+00	;
WV[15].omega[86] =	-3.473577E-03	;
WV[15].omega[87] =	1.482908E+00	;
WV[15].omega[88] =	1.991855E+00	;
WV[15].omega[89] =	-1.186099E+00	;
WV[15].omega[90] =	1.911730E+00	;
WV[15].w1[0][0] =	-1.728066E+01	;
WV[15].w1[0][1] =	2.210435E+00	;
WV[15].w1[0][2] =	2.348677E+00	;
WV[15].w1[0][3] =	2.812618E-01	;
WV[15].w1[0][4] =	2.417725E+00	;
WV[15].w1[0][5] =	1.246956E+00	;
WV[15].w1[0][6] =	-6.094876E-01	;
WV[15].w1[0][7] =	-5.309743E+00	;
WV[15].w1[0][8] =	-2.152382E+00	;
WV[15].w1[0][9] =	1.873751E+00	;
WV[15].w1[0][10] =	7.651818E-02	;
WV[15].w1[0][11] =	-1.529941E+01	;
WV[15].w1[0][12] =	1.839610E+00	;
WV[15].w1[0][13] =	-2.523251E+00	;
WV[15].w1[0][14] =	2.428332E+00	;
WV[15].w1[1][0] =	-2.273623E-02	;
WV[15].w1[1][1] =	-1.737884E+00	;
WV[15].w1[1][2] =	7.905896E-01	;
WV[15].w1[1][3] =	6.848911E-01	;
WV[15].w1[1][4] =	-3.599874E+00	;
WV[15].w1[1][5] =	-6.741247E-02	;
WV[15].w1[1][6] =	9.020591E-01	;
WV[15].w1[1][7] =	-7.354806E-03	;
WV[15].w1[1][8] =	9.320515E-02	;
WV[15].w1[1][9] =	1.300522E+00	;
WV[15].w1[1][10] =	-7.692179E-01	;
WV[15].w1[1][11] =	2.020262E+00	;
WV[15].w1[1][12] =	1.457440E+00	;
WV[15].w1[1][13] =	3.024908E+00	;
WV[15].w1[1][14] =	-4.385509E+00	;
WV[15].w1[2][0] =	5.227908E-01	;
WV[15].w1[2][1] =	6.498177E-01	;
WV[15].w1[2][2] =	-3.313994E-01	;
WV[15].w1[2][3] =	6.635243E-01	;
WV[15].w1[2][4] =	-3.528461E-01	;
WV[15].w1[2][5] =	-3.097976E-01	;
WV[15].w1[2][6] =	6.269263E-01	;
WV[15].w1[2][7] =	1.979896E-01	;
WV[15].w1[2][8] =	-5.230736E-01	;
WV[15].w1[2][9] =	-5.516489E-01	;
WV[15].w1[2][10] =	-6.330952E-01	;
WV[15].w1[2][11] =	1.594426E+01	;
WV[15].w1[2][12] =	-6.062846E-01	;
WV[15].w1[2][13] =	3.653671E-01	;
WV[15].w1[2][14] =	-3.582658E-01	;
WV[15].w1[3][0] =	2.421305E+00	;
WV[15].w1[3][1] =	7.068348E-01	;
WV[15].w1[3][2] =	-2.395018E+00	;
WV[15].w1[3][3] =	2.980795E+00	;
WV[15].w1[3][4] =	-1.885558E+00	;
WV[15].w1[3][5] =	-2.182956E+00	;
WV[15].w1[3][6] =	3.518668E+00	;
WV[15].w1[3][7] =	1.327941E+00	;
WV[15].w1[3][8] =	-2.535836E+00	;
WV[15].w1[3][9] =	-4.720176E+00	;
WV[15].w1[3][10] =	-3.051918E+00	;
WV[15].w1[3][11] =	-6.070693E+00	;
WV[15].w1[3][12] =	-5.186735E+00	;
WV[15].w1[3][13] =	2.165606E+00	;
WV[15].w1[3][14] =	-1.691555E+00	;
WV[15].w1[4][0] =	1.399478E+01	;
WV[15].w1[4][1] =	2.389718E+00	;
WV[15].w1[4][2] =	1.780490E+00	;
WV[15].w1[4][3] =	2.466100E+00	;
WV[15].w1[4][4] =	-2.008582E+00	;
WV[15].w1[4][5] =	2.876683E-01	;
WV[15].w1[4][6] =	3.439111E+00	;
WV[15].w1[4][7] =	3.260034E+00	;
WV[15].w1[4][8] =	-3.673472E+00	;
WV[15].w1[4][9] =	-1.727469E+00	;
WV[15].w1[4][10] =	-2.751030E+00	;
WV[15].w1[4][11] =	2.712121E+00	;
WV[15].w1[4][12] =	-1.946774E+00	;
WV[15].w1[4][13] =	2.165539E+00	;
WV[15].w1[4][14] =	-1.973277E+00	;
WV[15].w2[0][0] =	-1.152547E-01	;
WV[15].w2[1][0] =	-4.966613E-01	;
WV[15].w2[2][0] =	-5.353258E+00	;
WV[15].w2[3][0] =	1.614083E+00	;
WV[15].w2[4][0] =	3.080289E+00	;
WV[15].w2[5][0] =	3.921723E+00	;
WV[15].w2[6][0] =	1.062526E+00	;
WV[15].w2[7][0] =	-2.159306E-01	;
WV[15].w2[8][0] =	-8.224340E-01	;
WV[15].w2[9][0] =	-1.889153E+00	;
WV[15].w2[10][0] =	2.619577E+00	;
WV[15].w2[11][0] =	-3.473577E-03	;
WV[15].w2[12][0] =	1.482908E+00	;
WV[15].w2[13][0] =	1.991855E+00	;
WV[15].w2[14][0] =	-1.186099E+00	;
WV[15].w2[15][0] =	1.911730E+00	;
nu_star_MAX[15] =	1.00E+03	;
nu_star_MIN[15] =	1.00E-06	;
D_star_MAX[15]  =	1.00E+05	;
D_star_MIN[15]  =	1.00E-04	;
GGGMAX[15]      =	3.000000E-01	;
#endif

	return 0;
}
void SET_ARADI_IOTA(double RAX, double beta, double rho, double *paradi, double *piota)
{
	int iRAX;
	/* setteing of aradi and iota */
	iRAX=(int)(RAX * 100.0);
	switch(iRAX){
	case 345:if(beta==0.0){*paradi = 0.482714021199348;*piota = 0.52273+0.24083*pow(rho,2)+0.38*pow(rho,4)-0.30172*pow(rho,6)+0.32064*pow(rho,8);break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 350:if(beta==0.0){*paradi = 0.534982422591407;*piota = 0.45367+0.31539*pow(rho,2)+0.047693*pow(rho,4)+0.25885*pow(rho,6);break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 353:if(beta==0.0){*paradi = 0.512655337992553;*piota = 0.41827+0.45099*pow(rho,2)-0.13818*pow(rho,4)+0.56077*pow(rho,6);break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	case 360:/*if(beta==0.0){*paradi = 0.54;*piota = (double)0.33895+0.80324*pow(rho,2)-0.80193*pow(rho,4)+1.2353*pow(rho,6);break;}else*/
			 if(beta>=0.0&&beta<0.005){*paradi = 0.553836353538567;*piota = 0.36958+0.36808*pow(rho,2)+0.07362*pow(rho,4)+0.22762*pow(rho,6);break;}else
			 if(beta>=0.005&&beta<0.01){*paradi = 0.523247895;*piota = 0.35278+0.36803*pow(rho,2)+0.98597*pow(rho,4)-1.3195*pow(rho,6)+1.0837*pow(rho,8);break;}else 
			 if(beta>=0.01&&beta<0.02){*paradi = 0.53241591;*piota = 0.36715+0.16698*pow(rho,2)+1.2945*pow(rho,4)-1.4404*pow(rho,6)+0.96717*pow(rho,8);break;}else 
			 if(beta>=0.02&&beta<0.03){*paradi = 0.549259975;*piota = 0.51922-0.21291*pow(rho,2)+0.06784*pow(rho,4)+0.98405*pow(rho,6)-0.3254*pow(rho,8);break;}else 
			 if(beta=0.03){*paradi = 0.562574468180531;*piota = 0.72722-1.0365*pow(rho,2)+3.29*pow(rho,4)-6.721*pow(rho,6)+4.3923*pow(rho,8);break;}else 
			 {printf("(4)An improper value was substituted for beta:%f. \n", beta);exit(1);}
	case 375:if(beta>=0.00 && beta<0.01){*paradi = 0.563784329834497;*piota = 0.80631+0.050122*pow(rho,2)-1.9036*pow(rho,4)+1.3135*pow(rho,6);break;}else
			 if(beta>=0.01 && beta<0.02){*paradi = 0.537960220672511;*piota = 0.48537-0.10222*pow(rho,2)+0.66012*pow(rho,4)-0.34412*pow(rho,6)+0.36293*pow(rho,8)-0.078767*pow(rho,10);break;}else 
			 if(beta>=0.0199 && beta<0.03){*paradi = 0.540487355237937;*piota = 0.67798-0.56279*pow(rho,2)+1.7182*pow(rho,4)-4.075*pow(rho,6)+4.104*pow(rho,8)-1.159*pow(rho,10);break;}else 
			 if(beta>=0.0298 && beta<=0.05){*paradi = 0.540722689128731;*piota = 0.84479-0.42051*pow(rho,2)-1.9736*pow(rho,4)+8.1163*pow(rho,6)-13.335*pow(rho,8)+7.1142*pow(rho,10);break;}else 
			 {printf("An improper value was substituted for beta. \n");exit(1);}
	case 390:if(beta==0.0){*paradi = 0.560289757155521;*piota = 0.80266-0.093357*pow(rho,2)-0.93538*pow(rho,4)+0.73368*pow(rho,6);break;}else{printf("An improper value was substituted for beta. \n");exit(1);}
	default :printf("An improper value was substituted for RAX. \n");exit(1);
	}

	return;
}
