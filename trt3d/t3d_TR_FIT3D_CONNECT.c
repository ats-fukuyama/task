/*
 *  TR_FIT3D_CONNECT.c
 *  Y_task3D_modified
 *
 *  Created by WAKASA Arimitsu on 11/09/15.
 *  Copyright 2011 Arimitsu Wakasa MC. All rights reserved.
 *  最小二乗法：http://www.geocities.jp/supermisosan/saisyounizyouhou2.html
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <unistd.h>

#define NUM_COEF 7           //多項式の次数：(NUM_COEF-1)次
#define NUM_UNKNOWN_GAUSS NUM_COEF
#define NUM_DATA 60          // 多項式近似に用いるデータセット数

#define FIT3D_PATH "."
#define TASK3D_PATH  "."

#define BOZDATA_PATH "../profile_data/boz_datadir/newboz_a_109081t4240.dat"
#define BOZDATA_FORMAT "TXT"

#define TASK3D_PATH_LINK "ln_working_dir"
#define BOZDATA_PATH_LINK "ln_boz_data"

#define NBI_1 1 // ON=1 or OFF=0
#define NBI_2 1 // ON=1 or OFF=0
#define NBI_3 1 // ON=1 or OFF=0
#define NBI_4 1 // ON=1 or OFF=0
#define NBI_5 1 // ON=1 or OFF=0

#define NUM_NBI_UFILE 29

extern int task3d_tentative_flags_mp_ntpolynomial2ufile_flg_;
extern int task3d_tentative_flags_mp_fit3dtoufile_flg_;

extern void ftrd20_(void);
extern void depsum_(void);
extern void bozdata_translation_yt2mb_(char[],int);

void Gaussian_elimination_method(double [NUM_UNKNOWN_GAUSS][NUM_UNKNOWN_GAUSS+1],double [NUM_UNKNOWN_GAUSS]);
void polynomial_approximation(int,double [],double [], double[]);

void tr_fit3d_connect_(char* SHOT_NUM_tmp, double* time,int* nrmax, double* BB, double* rax,  double rg[], double rte[], double rti[], double rnn[],double pex_nbi_e[], double pex_nbi_i[]){
    
    int i, j,icnt,nr;
    char nbi_onoff[7]="000000\n";
    char bozdata_path2fort[1000]="";
    int len_bozdata_path2fort;
//    int len_SHOT_NUM, cnt_len_SHOT_NUM;
    char SHOT_NUM[80]="";
    double enbi1,enbi2,enbi3,enbi4,enbi5;
    double pnbi1,pnbi2,pnbi3,pnbi4,pnbi5;
    double SUM_qnbi_e[NUM_NBI_UFILE]={0.e0}, SUM_qnbi_i[NUM_NBI_UFILE]={0.e0};
    double qnbi_e[6][NUM_NBI_UFILE]={0.e0}, qnbi_i[6][NUM_NBI_UFILE]={0.e0},rg_nbi[NUM_NBI_UFILE];    
    double tmp;
    double coef_polynomial_PEXe[NUM_COEF],coef_polynomial_PEXi[NUM_COEF];
    double coef_polynomial_te[NUM_COEF]={
        2.4447,
        -2.6791,
        2.0942,
        -1.5424
    };
    double coef_polynomial_ti[NUM_COEF]={
        1.6922,
        -2.8675,
        2.9053,
        -1.4088
    };
    double coef_polynomial_nn[NUM_COEF]={
        1.6084,
        -1.0455,
        3.7143,
        -2.9011
    };
    
    /* Ajust the length of string (SHOT_NUM) */
    for(i=0;i<80;i++){
        if(*(SHOT_NUM_tmp+i)!=' '){
            SHOT_NUM[i]=*(SHOT_NUM_tmp+i);
        }
        else {
            SHOT_NUM[i]='\0';
            break;
        }
    }
    
/* Settel the NBI Power and particle energy  */    
/* Copy:below */
//    inbi1='on'
    enbi1=192.0996; pnbi1=6.0631;
//    inbi2='on'
    enbi2=177.5854; pnbi2=4.4140;
//    inbi3='on'
    enbi3=179.2896; pnbi3=4.7053;
//    inbi4='on'
    enbi4= 39.0014; pnbi4=5.3168;
//    inbi5='on'
    enbi5= 42.2265; pnbi5=5.7092;
    
/* Copy:above */    

    if(task3d_tentative_flags_mp_ntpolynomial2ufile_flg_==0){
        if(NUM_DATA != *nrmax){
            printf("\n\n");
            printf("===============================================================\n");
            printf("ERROR @ TR_FIT3D_CONNECT.c \n");
            printf("NUM_DATA %d != nrmax %d.\n",NUM_DATA,*nrmax);
            printf("Change the value of NUM_DATA (defined @ TR_FIT3D_CONNECT.c) =>> nrmax.\n\n");
            exit(0);
        }
        
        
        polynomial_approximation(*nrmax,rg,rte,coef_polynomial_te);
        polynomial_approximation(*nrmax,rg,rti,coef_polynomial_ti);
        polynomial_approximation(*nrmax,rg,rnn,coef_polynomial_nn);
    
        for(i=0;i<NUM_COEF;i++){
            coef_polynomial_te[i]=coef_polynomial_te[i]*1.e3;
            coef_polynomial_ti[i]=coef_polynomial_ti[i]*1.e3;
            coef_polynomial_nn[i]=coef_polynomial_nn[i]*1.e13; 
        }
    }
    else if(task3d_tentative_flags_mp_ntpolynomial2ufile_flg_==1 && task3d_tentative_flags_mp_fit3dtoufile_flg_ == 1){
        printf("\n");
        printf("===============================================================\n");
        printf("======   nT profiles are made by Polynomial eq. and    ========\n");
        printf("======   Translate to Ufile format.                    ========\n");
        printf("======     Te   : Handover_Te_Poly2UFILE.txt           ========\n");
        printf("======     Ti   : Handover_Ti_Poly2UFILE.txt           ========\n");
        printf("======   density: Handover_nn_Poly2UFILE.txt           ========\n");
        printf("===============================================================\n");
        printf("\n");
        for(i=3;i>=0;i--){
            coef_polynomial_te[i+i]=coef_polynomial_te[i]*1.e3;
            coef_polynomial_ti[i+i]=coef_polynomial_ti[i]*1.e3;
            coef_polynomial_nn[i+i]=coef_polynomial_nn[i]*1.e13; 
        }
        for(i=1;i<NUM_COEF;i=i+2){
            coef_polynomial_te[i]=0e0;
            coef_polynomial_ti[i]=0e0;
            coef_polynomial_nn[i]=0e0;             
        }
        for(nr=0;nr<*nrmax;nr++){
            rte[nr]=0.0;
            rti[nr]=0.0;
            rnn[nr]=0.0;
            for(i=0;i<NUM_UNKNOWN_GAUSS;i++) {
                rte[nr]+=coef_polynomial_te[i]*pow(rg[nr],i);
                rti[nr]+=coef_polynomial_ti[i]*pow(rg[nr],i);
                rnn[nr]+=coef_polynomial_nn[i]*pow(rg[nr],i);
            }
        }
        for(nr=0;nr<*nrmax;nr++){
            rte[nr]=rte[nr]/1e3;//unit>>[keV]
            rti[nr]=rti[nr]/1e3;//unit>>[keV]
            rnn[nr]=rnn[nr]/1e14;//unit>>[10^20m-3]
        }
    }
    else {
        printf("\n");
        printf("==================================================================\n");
        printf("======   ERROR @ FIT3D flags.                             ========\n");
        printf("======   Only if \"ntpolynomial2ufile_flg == 1\",           ========\n");
        printf("======   \"fit3dtoufile_flg == 1\" is allowed.              ========\n");
        printf("======   Check these flgs @ \"task3d_tentative_flags.f90\". ========\n");
        printf("==================================================================\n");
        printf("\n");
        exit(0);
    }


//#define CHCK_POLYNOMIAL_RNT
#ifdef CHCK_POLYNOMIAL_RNT
    double te_check[NUM_DATA],ti_check[NUM_DATA],nn_check[NUM_DATA];
    for(nr=0;nr<*nrmax;nr++){
        te_check[nr]=0.0;
        ti_check[nr]=0.0;
        nn_check[nr]=0.0;
        for(i=0;i<NUM_UNKNOWN_GAUSS;i++) {
            te_check[nr]+=coef_polynomial_te[i]*pow(rg[nr],i);
            ti_check[nr]+=coef_polynomial_ti[i]*pow(rg[nr],i);
            nn_check[nr]+=coef_polynomial_nn[i]*pow(rg[nr],i);
        }
    }
    FILE * RTN_CHECK;
    printf("\n");
    RTN_CHECK=fopen("RNT_CHECK.txt","w");
    fprintf(RTN_CHECK,"%4s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\t%20s\n","nr","time","rg[nr]","rte[nr]","rti[nr]","rnn[nr]","te_check[nr]","ti_check[nr]","nn_check[nr]");        
    for(nr=0;nr<*nrmax;nr++){
        fprintf(RTN_CHECK,"%4d\t%20.12e\t%20.12e\t%20.12e\t%20.12e\t%20.12e\t%20.12e\t%20.12e\t%20.12e\n",nr,*time,rg[nr],rte[nr],rti[nr],rnn[nr],te_check[nr]/1e3,ti_check[nr]/1e3,nn_check[nr]/1e19);        
        printf("%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",nr,*time,rg[nr],rte[nr],rti[nr],rnn[nr],te_check[nr]/1e3,ti_check[nr]/1e3,nn_check[nr]/1e19);        
    }
    fclose(RTN_CHECK);
    exit(0);
#endif
    
    printf("\n\nBOZ data check!\n");
    if(BOZDATA_FORMAT=="TXT"){
        printf("Translate from Yt to Mb.\n");
        strcpy(bozdata_path2fort,BOZDATA_PATH);
        len_bozdata_path2fort=strlen(bozdata_path2fort);
        bozdata_path2fort[len_bozdata_path2fort]=' ';
        printf("BOZ_data PATH: %s (%d characters)\n", bozdata_path2fort,len_bozdata_path2fort);
        bozdata_translation_yt2mb_(bozdata_path2fort, len_bozdata_path2fort);
    }
    else if(BOZDATA_FORMAT=="BIN"){
        printf("unneed Translation.\n");
    }
    else {
        printf("BOZ data file FORMAT ERROR!!\n");
        exit(1);
    }
    
//========================================================================================//
//========   Maiking Fit3d input files    ================================================//
//========================================================================================//
    
    char FIT_INP_FILENAME[150];
    char FIT_INP_FILENAMEwithPATH[200];
    char ECHO_FIT3D[1000];
    char ECHO_LN[1000];

    char char_time[6];
    FILE* FITINP;
    sprintf(FIT_INP_FILENAME, "%s%08.4f%s","t",*time,"ms.in");
    sprintf(FIT_INP_FILENAMEwithPATH, "%s%s%08.4f%s",TASK3D_PATH,"/t",*time,"ms.in");
    
    printf("Now making FIT3D input files @ %s.\n",FIT_INP_FILENAMEwithPATH);

    FITINP=fopen(FIT_INP_FILENAMEwithPATH,"w");
    sprintf(char_time,"%d%d%d%d%d",(int)(*time*1e0)%(int)10,(int)(*time*1e1)%(int)10,(int)(*time*1e2)%(int)10,(int)(*time*1e3)%(int)10,(int)(*time*1e4)%(int)10);

    
    printf("\nchar_time :: %s\n",char_time);
//    printf("%s\n",char_time);
    
    for(i=1;i<6;i++){
        if(NBI_1 == 1)nbi_onoff[1]='1';
        if(NBI_2 == 1)nbi_onoff[2]='1';
        if(NBI_3 == 1)nbi_onoff[3]='1';
        if(NBI_4 == 1)nbi_onoff[4]='1';
        if(NBI_5 == 1)nbi_onoff[5]='1';
    }
    
    fprintf(FITINP,"&cnfpls\n");
    fprintf(FITINP,"chshot= '%6s',\n",SHOT_NUM);
    fprintf(FITINP,"chtime= '%5s',\n",char_time);
    fprintf(FITINP,"raxcg = %9.3e,\n",*rax);
    fprintf(FITINP,"bax = %9.3e,\n",*BB);
    fprintf(FITINP,"ipg   = %c,\n",'0');
    fprintf(FITINP,"&end\n");
    fprintf(FITINP,"&cnfnbi\n");
    fprintf(FITINP,"inbi1=%c, enbi1=%10.4f, pnbi1=%10.4f\n",nbi_onoff[1],enbi1,pnbi1);
    fprintf(FITINP,"inbi2=%c, enbi2=%10.4f, pnbi2=%10.4f,\n",nbi_onoff[2],enbi2,pnbi2);
    fprintf(FITINP,"inbi3=%c, enbi3=%10.4f, pnbi3=%10.4f,\n",nbi_onoff[3],enbi3,pnbi3);
    fprintf(FITINP,"inbi4=%c, enbi4=%10.4f, pnbi4=%10.4f,\n",nbi_onoff[4],enbi4,pnbi4);
    fprintf(FITINP,"inbi5=%c, enbi5=%10.4f, pnbi5=%10.4f,\n" ,nbi_onoff[5],enbi5,pnbi5);
    fprintf(FITINP,"&end\n");
    fprintf(FITINP,"&prof\n");
    fprintf(FITINP,"anexm=\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\n",coef_polynomial_nn[0],coef_polynomial_nn[1],coef_polynomial_nn[2],coef_polynomial_nn[3],coef_polynomial_nn[4],coef_polynomial_nn[5],coef_polynomial_nn[6]);
    fprintf(FITINP,"atexm=\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\n",coef_polynomial_te[0],coef_polynomial_te[1],coef_polynomial_te[2],coef_polynomial_te[3],coef_polynomial_te[4],coef_polynomial_te[5],coef_polynomial_te[6]);
    fprintf(FITINP,"atixm=\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\t,\t%9.3E\n",coef_polynomial_ti[0],coef_polynomial_ti[1],coef_polynomial_ti[2],coef_polynomial_ti[3],coef_polynomial_ti[4],coef_polynomial_ti[5],coef_polynomial_ti[6]);
    fprintf(FITINP,"&end\n");
    fclose(FITINP);
    
    //========================================================================================//
    //========   Maiking Symbolic link                 =======================================//
    //========    1) TASK3D folder (working directry)  =======================================//
    //========    2) BOZ data folder                   =======================================//
    //========================================================================================//
    
    
    printf("\nTASK3D_PATH :: ");
    printf("%s\n",TASK3D_PATH);
    

    sprintf(ECHO_LN,"%s", "rm -fr ");
    strcat(ECHO_LN, TASK3D_PATH_LINK);
    printf("%s\n",ECHO_LN);
    system(ECHO_LN);    
    sprintf(ECHO_LN,"%s", "rm -fr ");
    strcat(ECHO_LN, BOZDATA_PATH_LINK);
    printf("%s\n",ECHO_LN);
    system(ECHO_LN);    

    sprintf(ECHO_LN,"%s","");
    strcat(ECHO_LN,"ln -s ");
    strcat(ECHO_LN,TASK3D_PATH);
    strcat(ECHO_LN," ");
    strcat(ECHO_LN,TASK3D_PATH_LINK);
    printf("%s\n",ECHO_LN);
    system(ECHO_LN);    
    sprintf(ECHO_LN,"%s","");
    strcat(ECHO_LN,"ln -s ");
    if(BOZDATA_FORMAT=="BIN"){
        strcat(ECHO_LN,BOZDATA_PATH);
    }
    else if(BOZDATA_FORMAT=="TXT"){
        strcat(ECHO_LN,TASK3D_PATH);
        strcat(ECHO_LN,"/BOZdata_mb_format.dat");
    }
    else {
        printf("BOZ data file LINK error. @ %s, %d\n",__FILE__,__LINE__);
        exit(1);
    }
    strcat(ECHO_LN," ");
    strcat(ECHO_LN,BOZDATA_PATH_LINK);
    printf("%s\n",ECHO_LN);
    system(ECHO_LN);    
    
    sprintf(ECHO_FIT3D,"%s", "echo \"");
    strcat(ECHO_FIT3D,FIT_INP_FILENAME);
    strcat(ECHO_FIT3D,"\n");
    strcat(ECHO_FIT3D,TASK3D_PATH_LINK);
    strcat(ECHO_FIT3D,"\n");
    strcat(ECHO_FIT3D,"1\n");
    strcat(ECHO_FIT3D,"1\n");
    strcat(ECHO_FIT3D,BOZDATA_PATH_LINK);
    strcat(ECHO_FIT3D,"\n\" | ./FITmake.sh");
    
    printf("START maing echo\n\n");
    printf("%s\n",ECHO_FIT3D);
    system(ECHO_FIT3D);
    printf("\nCompleted making FIT3D .sh files\n");
    
    char MCNBI_PATH[150];
    char ECHO_MCNBI[1000];    

    sprintf(MCNBI_PATH,"s%s@t%d%d%d%d",SHOT_NUM,(int)(*time*1e0)%(int)10,(int)(*time*1e1)%(int)10,(int)(*time*1e2)%(int)10,(int)(*time*1e3)%(int)10);

    //    printf("%s\n",MCNBI_PATH);
    if(NBI_1 == 1){
        sprintf(ECHO_MCNBI,"%s","./");
        strcat(ECHO_MCNBI,MCNBI_PATH);
        strcat(ECHO_MCNBI,"/sh/");
        strcat(ECHO_MCNBI,MCNBI_PATH);
        strcat(ECHO_MCNBI,"ctr1.sh");
        printf("NBI #1 \n");
        printf("%s\n",ECHO_MCNBI);
        system(ECHO_MCNBI);
    }
    if(NBI_2 == 1){
        sprintf(ECHO_MCNBI,"%s","./");
        strcat(ECHO_MCNBI,MCNBI_PATH);
        strcat(ECHO_MCNBI,"/sh/");
        strcat(ECHO_MCNBI,MCNBI_PATH);
        strcat(ECHO_MCNBI,"co2.sh");
        printf("NBI #2 \n");
        printf("%s\n",ECHO_MCNBI);
        system(ECHO_MCNBI);
    }
    if(NBI_3 == 1){
        sprintf(ECHO_MCNBI,"%s","./");
        strcat(ECHO_MCNBI,MCNBI_PATH);
        strcat(ECHO_MCNBI,"/sh/");
        strcat(ECHO_MCNBI,MCNBI_PATH);
        strcat(ECHO_MCNBI,"ctr3.sh");
        printf("NBI #3 \n");
        printf("%s\n",ECHO_MCNBI);
        system(ECHO_MCNBI);
    }
    if(NBI_4 == 1){
        sprintf(ECHO_MCNBI,"%s","./");
        strcat(ECHO_MCNBI,MCNBI_PATH);
        strcat(ECHO_MCNBI,"/sh/");
        strcat(ECHO_MCNBI,MCNBI_PATH);
        strcat(ECHO_MCNBI,"pb4.sh");
        printf("NBI #4 \n");
        printf("%s\n",ECHO_MCNBI);
        system(ECHO_MCNBI);
    }
    if(NBI_5 == 1){
        sprintf(ECHO_MCNBI,"%s","./");
        strcat(ECHO_MCNBI,MCNBI_PATH);
        strcat(ECHO_MCNBI,"/sh/");
        strcat(ECHO_MCNBI,MCNBI_PATH);
        strcat(ECHO_MCNBI,"pb5.sh");
        printf("NBI #5 \n");
        printf("%s\n",ECHO_MCNBI);
        system(ECHO_MCNBI);
    }
    
    printf("\nCompleted calculating MCNBI. \n");
    
    char FIT_FILENAME5[100];
    char FIT_FILENAME6[100];
    char FIT_FILENAME10[100];
    char FIT_FILENAME11[100];
    char FIT_FILENAME20[100];
    char FIT_FILENAME100[100];

    
    printf("NOW calculating NBI #1\n");
    if(NBI_1 == 1){
        sprintf(FIT_FILENAME5,"%s%s%s","fit_",MCNBI_PATH,"ctr1.in");
        sprintf(FIT_FILENAME6,"%s%s%s","fit_",MCNBI_PATH,"ctr1.out");
        sprintf(FIT_FILENAME10,"%s%s%s","fit_",MCNBI_PATH,"ctr1.out10");
        sprintf(FIT_FILENAME11,"%s%s%s","fit_",MCNBI_PATH,"ctr1.out11");
        sprintf(FIT_FILENAME20,"%s%s",MCNBI_PATH,"ctr1.out20");
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME20," ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME5," ","prof.in5");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        printf("NOW waiting 3sec for NBI#1 ftrd20();.\n");
        sleep(3);
        ftrd20_();
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out6"," ",MCNBI_PATH,"/",FIT_FILENAME6);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out10"," ",MCNBI_PATH,"/",FIT_FILENAME10);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out11"," ",MCNBI_PATH,"/",FIT_FILENAME11);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.in5");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
    } 

    printf("NOW calculating NBI #2\n");
    if(NBI_2 == 1){
        sprintf(FIT_FILENAME5,"%s%s%s","fit_",MCNBI_PATH,"co2.in");
        sprintf(FIT_FILENAME6,"%s%s%s","fit_",MCNBI_PATH,"co2.out");
        sprintf(FIT_FILENAME10,"%s%s%s","fit_",MCNBI_PATH,"co2.out10");
        sprintf(FIT_FILENAME11,"%s%s%s","fit_",MCNBI_PATH,"co2.out11");
        sprintf(FIT_FILENAME20,"%s%s",MCNBI_PATH,"co2.out20");
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME20," ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME5," ","prof.in5");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        printf("NOW waiting 3sec for NBI#2 ftrd20();.\n");
        sleep(3);
        ftrd20_();
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out6"," ",MCNBI_PATH,"/",FIT_FILENAME6);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out10"," ",MCNBI_PATH,"/",FIT_FILENAME10);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out11"," ",MCNBI_PATH,"/",FIT_FILENAME11);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.in5");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
    } 
    
    if(NBI_3 == 1){
        sprintf(FIT_FILENAME5,"%s%s%s","fit_",MCNBI_PATH,"ctr3.in");
        sprintf(FIT_FILENAME6,"%s%s%s","fit_",MCNBI_PATH,"ctr3.out");
        sprintf(FIT_FILENAME10,"%s%s%s","fit_",MCNBI_PATH,"ctr3.out10");
        sprintf(FIT_FILENAME11,"%s%s%s","fit_",MCNBI_PATH,"ctr3.out11");
        sprintf(FIT_FILENAME20,"%s%s",MCNBI_PATH,"ctr3.out20");
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME20," ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME5," ","prof.in5");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        printf("NOW waiting 3sec for NBI#3 ftrd20();.\n");
        sleep(3);        
        ftrd20_();
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out6"," ",MCNBI_PATH,"/",FIT_FILENAME6);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out10"," ",MCNBI_PATH,"/",FIT_FILENAME10);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out11"," ",MCNBI_PATH,"/",FIT_FILENAME11);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.in5");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
    } 

    if(NBI_4 == 1){
        sprintf(FIT_FILENAME5,"%s%s%s","fit_",MCNBI_PATH,"pb4.in");
        sprintf(FIT_FILENAME6,"%s%s%s","fit_",MCNBI_PATH,"pb4.out");
        
        sprintf(FIT_FILENAME10,"%s%s%s","fit_",MCNBI_PATH,"pb4.out10");
        sprintf(FIT_FILENAME11,"%s%s%s","fit_",MCNBI_PATH,"pb4.out11");
        sprintf(FIT_FILENAME20,"%s%s",MCNBI_PATH,"pb4.out20");
        sprintf(FIT_FILENAME100,"%s%s","fit_","pb4.out10");

        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME20," ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME5," ","prof.in5");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        printf("NOW waiting 3sec for NBI#4-1 ftrd20();.\n");
        sleep(3);
        
        ftrd20_();
        
        sprintf(ECHO_FIT3D,"%s%s%s%s","cp ","prof.out10"," ",FIT_FILENAME100);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out6"," ",MCNBI_PATH,"/",FIT_FILENAME6);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out10"," ",MCNBI_PATH,"/",FIT_FILENAME10);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out11"," ",MCNBI_PATH,"/",FIT_FILENAME11);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
//=========================
        sprintf(FIT_FILENAME10,"%s%s%s","fit_",MCNBI_PATH,"pb4.out20");
        sprintf(FIT_FILENAME11,"%s%s%s","fit_",MCNBI_PATH,"pb4.out21");
        sprintf(FIT_FILENAME20,"%s%s",MCNBI_PATH,"pb4.out30");
        sprintf(FIT_FILENAME100,"%s%s","fit_","pb4.out20");
        
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME20," ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        printf("NOW waiting 3sec for NBI#4-2 ftrd20();.\n");
        sleep(3);
        
        ftrd20_();
        
        sprintf(ECHO_FIT3D,"%s%s%s%s","cp ","prof.out10"," " ,FIT_FILENAME100);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);

        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out6"," ",MCNBI_PATH,"/",FIT_FILENAME6);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out10"," ",MCNBI_PATH,"/",FIT_FILENAME10);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out11"," ",MCNBI_PATH,"/",FIT_FILENAME11);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
//        =====================
        sprintf(FIT_FILENAME10,"%s%s%s","fit_",MCNBI_PATH,"pb4.out30");
        sprintf(FIT_FILENAME11,"%s%s%s","fit_",MCNBI_PATH,"pb4.out31");
        sprintf(FIT_FILENAME20,"%s%s",MCNBI_PATH,"pb4.out40");
        sprintf(FIT_FILENAME100,"%s%s","fit_","pb4.out30");
        
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME20," ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        printf("NOW waiting 3sec for NBI#4-3 ftrd20();.\n");
        sleep(3);
        
        ftrd20_();
        
        sprintf(ECHO_FIT3D,"%s%s%s%s","cp ","prof.out10"," ",FIT_FILENAME100);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);

        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out6"," ",MCNBI_PATH,"/",FIT_FILENAME6);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out10"," ",MCNBI_PATH,"/",FIT_FILENAME10);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out11"," ",MCNBI_PATH,"/",FIT_FILENAME11);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.in5");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);        
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);        
//      =====================
        printf("NOW waiting 3sec for NBI#4SUM depsum_();.\n");
        sleep(3);
        
//        depsum_();
        system("./depsum_100208.lm");
        
        sprintf(ECHO_FIT3D,"%s%s%s%s%s","mv fit_pb4.out40 ",MCNBI_PATH,"/fit_",MCNBI_PATH,"pb4.out40");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);

        sprintf(ECHO_FIT3D,"%s%s","rm ","fit_pb4.*");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);        
    } 
    
    if(NBI_5 == 1){
        sprintf(FIT_FILENAME5,"%s%s%s","fit_",MCNBI_PATH,"pb5.in");
        sprintf(FIT_FILENAME6,"%s%s%s","fit_",MCNBI_PATH,"pb5.out");
        
        sprintf(FIT_FILENAME10,"%s%s%s","fit_",MCNBI_PATH,"pb5.out10");
        sprintf(FIT_FILENAME11,"%s%s%s","fit_",MCNBI_PATH,"pb5.out11");
        sprintf(FIT_FILENAME20,"%s%s",MCNBI_PATH,"pb5.out20");
        sprintf(FIT_FILENAME100,"%s%s","fit_","pb4.out10");
        
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME20," ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME5," ","prof.in5");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        printf("NOW waiting 3sec for NBI#5-1 ftrd20();.\n");
        sleep(3);        
        ftrd20_();
        
        sprintf(ECHO_FIT3D,"%s%s%s%s","cp ","prof.out10"," ",FIT_FILENAME100);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out6"," ",MCNBI_PATH,"/",FIT_FILENAME6);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out10"," ",MCNBI_PATH,"/",FIT_FILENAME10);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out11"," ",MCNBI_PATH,"/",FIT_FILENAME11);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        //=========================
        sprintf(FIT_FILENAME10,"%s%s%s","fit_",MCNBI_PATH,"pb5.out20");
        sprintf(FIT_FILENAME11,"%s%s%s","fit_",MCNBI_PATH,"pb5.out21");
        sprintf(FIT_FILENAME20,"%s%s",MCNBI_PATH,"pb5.out30");
        sprintf(FIT_FILENAME100,"%s%s","fit_","pb4.out20");
        
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME20," ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        printf("NOW waiting 3sec for NBI#5-2 ftrd20();.\n");
        sleep(3);        
        
        ftrd20_();
        
        sprintf(ECHO_FIT3D,"%s%s%s%s","cp ","prof.out10"," " ,FIT_FILENAME100);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out6"," ",MCNBI_PATH,"/",FIT_FILENAME6);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out10"," ",MCNBI_PATH,"/",FIT_FILENAME10);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out11"," ",MCNBI_PATH,"/",FIT_FILENAME11);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        //        =====================
        sprintf(FIT_FILENAME10,"%s%s%s","fit_",MCNBI_PATH,"pb5.out30");
        sprintf(FIT_FILENAME11,"%s%s%s","fit_",MCNBI_PATH,"pb5.out31");
        sprintf(FIT_FILENAME20,"%s%s",MCNBI_PATH,"pb5.out40");
        sprintf(FIT_FILENAME100,"%s%s","fit_","pb4.out30");
        
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","cp ",MCNBI_PATH,"/",FIT_FILENAME20," ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        printf("NOW waiting 3sec for NBI#5-3 ftrd20();.\n");
        sleep(3);        
        
        ftrd20_();
        
        sprintf(ECHO_FIT3D,"%s%s%s%s","cp ","prof.out10"," ",FIT_FILENAME100);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out6"," ",MCNBI_PATH,"/",FIT_FILENAME6);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out10"," ",MCNBI_PATH,"/",FIT_FILENAME10);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s%s%s%s%s","mv ","prof.out11"," ",MCNBI_PATH,"/",FIT_FILENAME11);
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.in5");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);        
        sprintf(ECHO_FIT3D,"%s%s","rm ","prof.out20");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);        
        //      =====================
        printf("NOW waiting 3sec for NBI#5SUM depsum_();.\n");
        sleep(3);        
        
//        depsum_();
        system("./depsum_100208.lm");

        sprintf(ECHO_FIT3D,"%s%s%s%s%s","mv fit_pb4.out40 ",MCNBI_PATH,"/fit_",MCNBI_PATH,"pb5.out40");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);
        
        sprintf(ECHO_FIT3D,"%s%s","rm ","fit_pb4.*");
        printf("%s\n",ECHO_FIT3D);
        system(ECHO_FIT3D);        
    } 
    
    printf("\nCompleted execution of FIT3D.\n");
    printf("And next, read and sum up the injection Power, EX.\n\n");

    
    // ====================  READ from FIT3D output files =======================//
    
    char QNBI_FILENAME[150];
    FILE* FITOUT; 
     if(NBI_1 == 1){
         sprintf(QNBI_FILENAME,"%s/fit_%sctr1.out10",MCNBI_PATH,MCNBI_PATH);
         FITOUT=fopen(QNBI_FILENAME,"r");
         for(i=0;i<NUM_NBI_UFILE;i++){
             fscanf(FITOUT,"%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n",
                    &icnt,&rg_nbi[i],&tmp,&qnbi_e[1][i],&qnbi_i[1][i],
                    &tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp);
//             printf("%d\t%e\t%e\t%e\n",i,rg_nbi[i],qnbi_e[1][i],qnbi_i[1][i]);
         }
         fclose(FITOUT);
     }
     if(NBI_2 == 1){
         sprintf(QNBI_FILENAME,"%s/fit_%sco2.out10",MCNBI_PATH,MCNBI_PATH);
         FITOUT=fopen(QNBI_FILENAME,"r");
         for(i=0;i<NUM_NBI_UFILE;i++){
             fscanf(FITOUT,"%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n",
                    &icnt,&rg_nbi[i],&tmp,&qnbi_e[2][i],&qnbi_i[2][i],
                    &tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp);
             //             printf("%d\t%e\t%e\t%e\n",i,rg_nbi[i],qnbi_e[1][i],qnbi_i[1][i]);
         }
         fclose(FITOUT);
     }
     if(NBI_3 == 1){
         sprintf(QNBI_FILENAME,"%s/fit_%sctr3.out10",MCNBI_PATH,MCNBI_PATH);
         FITOUT=fopen(QNBI_FILENAME,"r");
         for(i=0;i<NUM_NBI_UFILE;i++){
             fscanf(FITOUT,"%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n",
                    &icnt,&rg_nbi[i],&tmp,&qnbi_e[3][i],&qnbi_i[3][i],
                    &tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp);
             //             printf("%d\t%e\t%e\t%e\n",i,rg_nbi[i],qnbi_e[1][i],qnbi_i[1][i]);
         }
         fclose(FITOUT);
     }
     if(NBI_4 == 1){
         sprintf(QNBI_FILENAME,"%s/fit_%spb4.out40",MCNBI_PATH,MCNBI_PATH);
         FITOUT=fopen(QNBI_FILENAME,"r");
         for(i=0;i<NUM_NBI_UFILE;i++){
             fscanf(FITOUT,"%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n",
                    &icnt,&rg_nbi[i],&tmp,&qnbi_e[4][i],&qnbi_i[4][i],
                    &tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp);
             //             printf("%d\t%e\t%e\t%e\n",i,rg_nbi[i],qnbi_e[1][i],qnbi_i[1][i]);
         }
         fclose(FITOUT);
     }
     if(NBI_5 == 1){
         sprintf(QNBI_FILENAME,"%s/fit_%spb5.out40",MCNBI_PATH,MCNBI_PATH);
         FITOUT=fopen(QNBI_FILENAME,"r");
         for(i=0;i<NUM_NBI_UFILE;i++){
             fscanf(FITOUT,"%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n",
                    &icnt,&rg_nbi[i],&tmp,&qnbi_e[5][i],&qnbi_i[5][i],
                    &tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp);
             //             printf("%d\t%e\t%e\t%e\n",i,rg_nbi[i],qnbi_e[1][i],qnbi_i[1][i]);
         }
         fclose(FITOUT);         
     }

    for(j=1;j<NUM_NBI_UFILE;j++){
        for(i=1;i<6;i++){
            SUM_qnbi_e[j]=SUM_qnbi_e[j]+qnbi_e[i][j];
            SUM_qnbi_i[j]=SUM_qnbi_i[j]+qnbi_i[i][j];
        }
    }
    
    polynomial_approximation(NUM_NBI_UFILE,rg_nbi,SUM_qnbi_e,coef_polynomial_PEXe);
    polynomial_approximation(NUM_NBI_UFILE,rg_nbi,SUM_qnbi_i,coef_polynomial_PEXi);
    for(nr=0;nr<*nrmax;nr++){
        pex_nbi_e[nr]=0.0;
        pex_nbi_i[nr]=0.0;
        for(i=0;i<NUM_UNKNOWN_GAUSS;i++) {
            pex_nbi_e[nr]+=coef_polynomial_PEXe[i]*pow(rg[nr],i)*1e6;
            pex_nbi_i[nr]+=coef_polynomial_PEXi[i]*pow(rg[nr],i)*1e6;
        }
    }
    
    FILE* FITCHECK;
    static int flg_FITCHECK=0;
    if(flg_FITCHECK == 0){
        FITCHECK=fopen("FIT3D_results_for_PEX.txt","w");
    }
    else {
        FITCHECK=fopen("FIT3D_results_for_PEX.txt","a");
        flg_FITCHECK=1;
    }
    fprintf(FITCHECK, "#i\trg_nbi[i]\tSUM_qnbi_e[i]\tqnbi_e[1][i]\tqnbi_e[2][i]\tqnbi_e[3][i]\tqnbi_e[4][i]\tqnbi_e[5][i]\tSUM_qnbi_i[i]\tqnbi_i[1][i]\tqnbi_i[2][i]\tqnbi_i[3][i]\tqnbi_i[4][i]\tqnbi_i[5][i]\n");
    for(i=0;i<NUM_NBI_UFILE;i++){
        fprintf(FITCHECK, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",i,rg_nbi[i],SUM_qnbi_e[i],qnbi_e[1][i],qnbi_e[2][i],qnbi_e[3][i],qnbi_e[4][i],qnbi_e[5][i],SUM_qnbi_i[i],qnbi_i[1][i],qnbi_i[2][i],qnbi_i[3][i],qnbi_i[4][i],qnbi_i[5][i]);
    }
    fprintf(FITCHECK, "\n");
    fclose(FITCHECK);
    
    
#define CHCK_POLYNOMIAL_RNT_PEX    
#ifdef CHCK_POLYNOMIAL_RNT_PEX
    for(nr=0;nr<NUM_NBI_UFILE;nr++){
        printf("%d\t%e\t%e\t%e\n",nr,rg_nbi[nr],SUM_qnbi_e[nr],SUM_qnbi_i[nr]);
    }
    printf("");
    
    printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\n",coef_polynomial_PEXe[0],coef_polynomial_PEXe[1],coef_polynomial_PEXe[2],coef_polynomial_PEXe[3],coef_polynomial_PEXe[4],coef_polynomial_PEXe[5],coef_polynomial_PEXe[6]);
    printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\n",coef_polynomial_PEXi[0],coef_polynomial_PEXi[1],coef_polynomial_PEXi[2],coef_polynomial_PEXi[3],coef_polynomial_PEXi[4],coef_polynomial_PEXi[5],coef_polynomial_PEXi[6]);
    printf("===================================\n");
    
    FILE * PEX_CHECK;
    PEX_CHECK=fopen("PEX_CHECK.txt","w");
    fprintf(PEX_CHECK,"%4s\t%20s\t%20s\t%20s\t%20s\n","nr","*time","rg[nr]","pex_nbi_e[nr]","pex_nbi_i[nr]");
    for(nr=0;nr<*nrmax;nr++){
        fprintf(PEX_CHECK,"%4d\t%20.12e\t%20.12e\t%20.12e\t%20.12e\n",nr,*time,rg[nr],pex_nbi_e[nr],pex_nbi_i[nr]);
    }
    fprintf(PEX_CHECK,"");
#endif
    
    return ;
}


void polynomial_approximation(int num_maxdata,double xx_basis[NUM_DATA],double yy_basis[NUM_DATA],double coef_polynomial[NUM_COEF])
{
    int i,j,k;
    double matrix_for_gauss[NUM_COEF][NUM_COEF+1];
    
    for(i=0;i<NUM_COEF;i++) {
        for(j=0;j<NUM_COEF+1;j++) {
            matrix_for_gauss[i][j]=0.0;
        }
    }
    
    for(i=0;i<NUM_COEF;i++) {
        for(j=0;j<NUM_COEF;j++) {
            for(k=0;k<num_maxdata;k++) {
                matrix_for_gauss[i][j]+=pow(xx_basis[k],i+j);
            }
        }
    }
    for(i=0;i<NUM_COEF;i++) {
        for(k=0;k<num_maxdata;k++) {
            matrix_for_gauss[i][NUM_COEF]+=pow(xx_basis[k],i)*yy_basis[k];
        }
    }
    Gaussian_elimination_method(matrix_for_gauss,coef_polynomial);
    
}
void Gaussian_elimination_method(double a[NUM_UNKNOWN_GAUSS][NUM_UNKNOWN_GAUSS+1],double xx[NUM_UNKNOWN_GAUSS])
{
    int i,j,k,l,pivot;
    double x[NUM_UNKNOWN_GAUSS];
    double p,q,m,b[1][NUM_UNKNOWN_GAUSS+1];
    
    for(i=0;i<NUM_UNKNOWN_GAUSS;i++) {
        m=0;
        pivot=i;
        
        for(l=i;l<NUM_UNKNOWN_GAUSS;l++) {
            if(fabs(a[l][i])>m) { 
                m=fabs(a[l][i]);
                pivot=l;
            }
        }
        
        if(pivot!=i) {                          
            for(j=0;j<NUM_UNKNOWN_GAUSS+1;j++) {
                b[0][j]=a[i][j];        
                a[i][j]=a[pivot][j];
                a[pivot][j]=b[0][j];
            }
        }
    }
    
    for(k=0;k<NUM_UNKNOWN_GAUSS;k++) {
        p=a[k][k];
        a[k][k]=1;
        
        for(j=k+1;j<NUM_UNKNOWN_GAUSS+1;j++) {
            a[k][j]/=p;
        }
        
        for(i=k+1;i<NUM_UNKNOWN_GAUSS;i++) {
            q=a[i][k];
            
            for(j=k+1;j<NUM_UNKNOWN_GAUSS+1;j++) {
                a[i][j]-=q*a[k][j];
            }
            a[i][k]=0;
        }
    }
    
    for(i=NUM_UNKNOWN_GAUSS-1;i>=0;i--) {
        x[i]=a[i][NUM_UNKNOWN_GAUSS];
        for(j=NUM_UNKNOWN_GAUSS-1;j>i;j--) {
            x[i]-=a[i][j]*x[j];
        }
    }
    for(i=0;i<NUM_COEF;i++) {
        xx[i]=x[i];
    }    
}

