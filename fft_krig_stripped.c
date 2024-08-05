#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
//#include <complex.h>
#include <complex.h>
#include <time.h>
#include <assert.h>
#include <libgen.h> // for basename()
#include "cspline.h"
#include <getopt.h>
#include "slalib.h"
#include "uvfits.h"
#include "primary_beam.h"
#include "fitsio.h"
#include <omp.h>
	//#include <gsl/gsl_blas.h>
#include <cblas.h>

// gcc -Wall -g lssa_fg.c cspline.c uvfits.c -o lssa_fg -L/usr/local/lib -llapack -lcblas -lm -lcfitsio

#define expand_factor 1
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

/* function prototypes for LAPACK */
int cpotrf_(char *uplo, int *n, float complex *a, int *lda, int *info) ;
int cpotri_(char *uplo, int *n, float complex *a, int *lda, int *info ) ;
int cpotrs_(char *uplo, int *n, int *nrhs, float complex *a, int *lda, float complex *b, int *ldb, int *info );
//int zgetrf_(int *n, int *m,  double complex *a, int *lda, int *ipiv, int *info) ;
int zgetrf_(int *n, int *m,  double complex *a, int *lda, int *ipiv, int *info) ;
int zgetrs_(char *trans, int *n, int *nrhs, double complex *a, int *lda, int *ipiv, double complex *b, int *ldb, int *info);
int zgetri_(int *n, double complex *a, int *lda, int *ipiv, double complex *work, int *lwork, int *info);

/* Function prototypes */
double bessj1( double x );
int lssa_cross_pro(float u, int Nchan, float complex* vis0, float complex* vis1, float complex* flags0, float complex* flags1, float* weights0, float* totpower, float* crosspower, float* flagpower, float* residpower, float* residpowerimag, float* crossvar, float sigma2, int bias_mode, int band, float t_sample);
int dft_cross_pro(float u, int Nchan, float complex* vis0, float complex* vis1, float complex* flags0, float complex* flags1, float* weights0, float* totpower, float* crosspower, float* flagpower, float* residpower, float* residpowerimag, float* crossvar, float sigma2, int bias_mode, int band, float t_sample);
void interlin(float complex* input, int sizein, float* input_loc, float* output_loc, int sizeout, float complex* output);
void interlinf(float* input, int sizein, float* input_loc, float* output_loc, int sizeout, float* output);
void interpad(float complex* input, int sizein, float* input_loc, float* output_loc, int sizeout, float complex* output);
void interpadf(float* input, int sizein, float* input_loc, float* output_loc, int sizeout, float* output);
void free_mem3D(double ***matrix, int xsize, int ysize, int zsize);
int omp_get_num_threads(void);
int minimum_location(float* array, int size);
double vectorVectorMultiply(int vsize, double vec1[], double vec2[]);
void free_memFloat(float **matrix, int xsize, int ysize);
int matrixVectorMultiply(int vsize, double vec[], double **mat, double outVector[]);
double **create2DMatrix(int sizex, int sizey);
float complex **create2DMatrixComplex(int sizex, int sizey);
double ***create3DMatrix(int sizex, int sizey, int sizez);
//float ***create3DMatrixFloat(int sizex, int sizey, int sizez);
float complex ***create3DMatrixComplex(int sizex, int sizey, int sizez);
int ***create3DMatrixInt(int sizex, int sizey, int sizez);
long ***create3DMatrixLong(int sizex, int sizey, int sizez);
double ****create4DMatrix(int sizex, int sizey, int sizez, int sizew);
void memzero2DFloat(float **matrix,int xsize, int ysize);
long ****create4DMatrixLong(int sizex, int sizey, int sizez, int sizew);
float **create2DMatrixFloat(int sizex, int sizey);
unsigned long **create2DMatrixShort(int sizex, int sizey);
int **create2DMatrixInt(int sizex, int sizey);
long int **create2DMatrixLong(int sizex, int sizey);
int matrixMatrixMultiply(int x1size, int y1size, int x2size, int y2size, double **mat1, double **mat2, double **outMat);
int CmatrixMatrixMultiply(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag);
int CmatrixMatrixMultiplyConj(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag);
int CmatrixMatrixMultiplyConjTrans(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag);
int matrixTranspose(int xsize, int ysize, double **mat_in, double **mat_out);
double *create1DVector(int size);
void where(long *vec_in, long int size, int limit, long int *vec_out, int *count);
double trace_mat(double **matrix, int size);
void memzero2D(double **matrix,int xsize, int ysize);
void memzero2Dint(int **matrix,int xsize, int ysize);
void memzero2Dlong(long int **matrix,int xsize, int ysize);
void memzero3D(double ***matrix,int xsize, int ysize, int zsize);
void memzero3Dint(int ***matrix,int xsize, int ysize, int zsize);
void free_mem(double **matrix, int xsize, int ysize);
void clean_up();

/* global variables */
float global_period = PERIOD;
float global_chanwidth = CHAN_WIDTH;

int debug=0;
char *infilename=NULL;
int flag_dc=0;
FILE *fpd;

void usage() {
    fprintf(stderr,"Usage:\n");
    fprintf(stderr,"test_readuvfits <options> -i inputfilename\n");
    fprintf(stderr,"\t-d debug_level (0=no debugging)\n");
    exit(0);
}

void parse_cmdline(const int argc,char * const argv[]) {
    int result=0;
    const char optstring[] = "d:i:p:c:f:";

    while ( (result = getopt(argc, argv, optstring)) != -1 ) {
        switch (result) {
          case 'd': debug = atoi(optarg);
            break;
          case 'i': infilename = optarg;
            break;
            case 'p': global_period = atof(optarg);
                break;
            case 'c': global_chanwidth = atof(optarg);
                break;
	case 'f': flag_dc = atoi(optarg);
                break;
          default:
              fprintf(stderr,"unknown option: %c\n",result);
              usage(argv[0]);
        }
    }
}





/* Main */
int main (int argc, char *argv[]) {

  int i,j,k,size_fin,Nchan, nbins, **fg_num=NULL, band=0,Nchanall,Nfreq,mm,Nst=19;
  FILE *fptr=NULL,*fptrv=NULL,*fpbb=NULL,*fpflag=NULL,*flog=NULL,*fptrg=NULL;
  float *uu=NULL, *vv=NULL, *lperp=NULL, ***weights0=NULL, ***weights1=NULL,**power=NULL,**powertemp=NULL,**powerflag=NULL,**powertot=NULL,**power_resid=NULL,**power_residimag=NULL,**lssa_weights=NULL,maxu=0.,***lssa_weights_cube=NULL,***powercube=NULL,***residcube=NULL,obsvolume=0.;
  double u_size,corr_length_factor=0.;
  float complex ***flags0=NULL, ***flags1=NULL, ***vis0=NULL, ***vis1=NULL,sigma2=0,tsample=16.,Tsys=440.,***tempvis1=NULL,***tempvis0=NULL,***tempflags0=NULL, ***tempflags1=NULL;
  float lmax, space=0.,**krig_weights=NULL,temp=0.,temp2=0.,***tempweights0=NULL, ***tempweights1=NULL;
  char *date=NULL,syslogfile[1024],*pol=NULL,*extension=NULL,bvfilename[1024],bbfilename[1024],filename_flag[1024],flagfilename[1024],outfilename1[1024],outfilename2[1024],outfilename3[1024],outfilename4[1024],outfilename5[1024],outfilename6[1024],outfilename1_1[1024],outfilename5_1[1024],lperpfilename[1024],covfilename[1024],covgsname[1024],outfilename5_2[1024],outfilename5_3[1024],outfilenameresid[1024],paramsfilename[2014];
  unsigned long dim1;
  int nthreads=1,bias_mode=0;
	
  nthreads = omp_get_max_threads();
  printf("max threads: %d\n",nthreads);

 //   nthreads=1;

  fpd = stderr;

// Information


if (argc < 2){
    fprintf(stderr,"Usage: <options> %s input_extension Nchan nbins pol maxu output_extension bias_mode (0/10/11/12/13) band \n",argv[0]);
    fprintf(stderr,"\t -p period\n");
    fprintf(stderr,"\t -c chanwidth\n");
    fprintf(stderr,"\t -f flag_dc\n");
    exit(1);
}


printf("Program to compute the LS spectral power using diff, tot and weights uvf binary files.\n\n");

date = argv[6];
Nchan = atoi(argv[2]);
nbins = atoi(argv[3]);
pol = argv[4];
maxu = atof(argv[5]);
extension = argv[1];
bias_mode = atoi(argv[7]);
band = atoi(argv[8]);
    
    parse_cmdline(argc,argv);

    
tsample = global_period;

Nfreq = Nchan;

/* Open log file */
	sprintf(syslogfile,"%ssyslog_lssa.txt",getenv("OUTPUTDIR"));
	printf("syslogfile: %s\n",syslogfile);
	if ((flog=fopen(syslogfile,"w")) == NULL){
        fprintf(stderr,"Cannot open output log file\n");
	return 1;
	}

    /* define uv grid size and spacing */

u_size = floor(maxu/DELTA_U);
    printf("u_size %g u_size+1 %g 2.*u_size %g\n",u_size,u_size+1,2.*u_size);
    

    
/* Setup lm bins --- linear bins */

 uu = calloc(u_size+1,sizeof(float));
 vv = calloc(2.*u_size,sizeof(float));
 for (i=0;i<u_size+1;i++) uu[i] = (i)/((double)u_size)*(maxu);
 for (i=0;i<2.*u_size;i++) vv[i] = (i-u_size)/((double)u_size)*(maxu);

space = (vv[1]-vv[0]);

lmax = maxu*1.1;
lperp = calloc(nbins+1, sizeof(float));
for (i=0;i<nbins+1;i++) lperp[i] = (float)i/((float)nbins+1.)*lmax;
    
 

	/* read-in input files */

	fprintf(flog,"Opening and reading diff files...\n");

	sprintf(bvfilename,"%svis_diff_%s.%s.dat",getenv("INPUTDIR"),pol,extension);
	sprintf(flagfilename,"%snoise_diff_%s.%s.dat",getenv("INPUTDIR"),pol,extension);
	sprintf(bbfilename,"%sweights_%s.%s.dat",getenv("INPUTDIR"),pol,extension);

	printf("bvfilename: %s\n",bvfilename);

	vis0 = create3DMatrixComplex((int)u_size+1,(int)2.*u_size,Nchan);
	flags0 = create3DMatrixComplex((int)u_size+1,(int)2.*u_size,Nchan);
	weights0 = create3DMatrixFloat((int)u_size+1,(int)2.*u_size,Nchan);
    
//    printf("Defined arrays\n");

	if ((fptrv = fopen(bvfilename,"r")) == NULL){
		// file does not exist - exit
		fprintf(flog,"Vis diff file does not exist... exiting...\n");
		return 1;
	}


    if ((fpflag = fopen(flagfilename,"r")) == NULL){
        // file does not exist - exit
        fprintf(flog,"Noise diff file does not exist... exiting...\n");
        return 1;
    }
    
printf("bvfilename: %s\n",bvfilename);
	if ((fpbb = fopen(bbfilename,"r")) == NULL){
		// file does not exist - exit
		fprintf(flog,"Weights file does not exist... exiting...\n");
		return 1;
	}

//    printf("bvfilename: %s\n",bvfilename);
    
	for (i=0;i<(int)u_size+1;i++){
 //       printf("i %d\n",i);
		for (j=0;j<(int)2.*u_size;j++){
		fread(vis0[i][j],sizeof(float complex),Nchan,fptrv);
		fread(flags0[i][j],sizeof(float complex),Nchan,fpflag);
		fread(weights0[i][j],sizeof(float),Nchan,fpbb);
		}
	}

	fclose(fptrv);
	fclose(fpflag);
	fclose(fpbb);

	printf("weights0: %f %f %f %f\n",weights0[2][100][0],weights0[4][100][0],weights0[2][130][0],weights0[23][140][0]);

	fprintf(flog,"Diff/weights files read...\n");

	/* Set 1 */

	fprintf(flog,"Opening and reading set 1 files...\n");

	sprintf(bvfilename,"%svis_tot_%s.%s.dat",getenv("INPUTDIR"),pol,extension);
	sprintf(flagfilename,"%snoise_tot_%s.%s.dat",getenv("INPUTDIR"),pol,extension);

	vis1 = create3DMatrixComplex(u_size+1,2.*u_size,Nchan);
	flags1 = create3DMatrixComplex(u_size+1,2.*u_size,Nchan);

	if ((fptrv = fopen(bvfilename,"r")) == NULL){
		// file does not exist - exit
		fprintf(flog,"Vis Tot file does not exist... exiting...\n");
		return 1;
	}
	if ((fpflag = fopen(flagfilename,"r")) == NULL){
		// file does not exist - exit
		fprintf(flog,"Noise Tot file does not exist... exiting...\n");
		return 1;
	}

	for (i=0;i<u_size+1;i++){
		for (j=0;j<2.*u_size;j++){
		fread(vis1[i][j],sizeof(float complex),Nchan,fptrv);
		fread(flags1[i][j],sizeof(float complex),Nchan,fpflag);
		}
	}


	fclose(fptrv);
	fclose(fpflag);

//    flags1 = vis1;
//   flags0 = vis0;
    
	fprintf(flog,"Tot files read...\n");
    


    
    printf("vis %g %g %g %g weights %f %f %f %f\n",creal(vis1[10][610][0]),creal(vis1[10][610][144]),creal(vis1[10][610][16]),creal(vis1[10][610][15]),weights0[10][610][0],weights0[10][610][144],weights0[10][610][16],weights0[10][610][15]);
    
  

    if ((bias_mode == 10)||(bias_mode == 13)||(bias_mode == 14)||(bias_mode == 15)||(bias_mode == 11)||(bias_mode == 12)||(bias_mode == 20)||(bias_mode == 21)||(bias_mode == 22)){
        
            //       Nchan = 232;
            //	Nst = 19;
        Nchan= 192;
        if (bias_mode == 10) Nst = 192;
        if (bias_mode == 11) Nst = 0;
        if (bias_mode == 12) Nst = 96;
        if (bias_mode == 13) Nst = 48;
        if (bias_mode == 14){
            Nst = 196;
            Nchan = 184;
        }
        if (bias_mode == 15){
            Nst = 52;
            Nchan=184;
        }
	if (bias_mode == 20) Nst = 192;
	if (bias_mode == 21) Nst = 0;
	if (bias_mode == 22) Nst = 96;
        Nfreq = Nchan;
        
        printf("Nchan: %d %d\n",Nchan,Nfreq);
        
        tempvis0 = create3DMatrixComplex((int)u_size+1,(int)2.*u_size,Nchan);
        tempflags0 = create3DMatrixComplex((int)u_size+1,(int)2.*u_size,Nchan);
        tempweights0 = create3DMatrixFloat((int)u_size+1,(int)2.*u_size,Nchan);
        
        for (i=0;i<u_size+1;i++){
            for (j=0;j<2.*u_size;j++){
                for (k=0;k<Nchan;k++){
         
                    tempvis0[i][j][k] = vis0[i][j][k+Nst];
       //             tempvis0[i][j][k] = 0.;
                    tempweights0[i][j][k] = weights0[i][j][k+Nst];
                    tempflags0[i][j][k] = flags0[i][j][k+Nst];
                    
                
                }
            }
        }
        
        vis0 = tempvis0;
        weights0 = tempweights0;
        flags0 = tempflags0;

        
        tempvis1 = create3DMatrixComplex((int)u_size+1,(int)2.*u_size,Nchan);
        tempflags1 = create3DMatrixComplex((int)u_size+1,(int)2.*u_size,Nchan);
//        tempweights1 = create3DMatrixFloat((int)u_size+1,(int)2.*u_size,Nchan);
        
        for (i=0;i<u_size+1;i++){
            for (j=0;j<2.*u_size;j++){
                for (k=0;k<Nchan;k++){
                    
                    tempvis1[i][j][k] = vis1[i][j][k+Nst];
            //        tempvis1[i][j][k] = 1.;
 //                   tempweights1[i][j][k] = weights1[i][j][k+19];
                    tempflags1[i][j][k] = flags1[i][j][k+Nst];
                    
                    
                }
            }
        }
        
        vis1 = tempvis1;
//        weights1 = tempweights1;
        flags1 = tempflags1;
        
        
    }
    
    printf("vis %g %g %g %g weights %f %f %f %f\n",creal(vis1[10][610][0]),creal(vis1[10][610][144]),creal(vis1[10][610][16]),creal(vis1[10][610][15]),weights0[10][610][0],weights0[10][610][144],weights0[10][610][16],weights0[10][610][15]);



	/* Allocate output vectors */

    int sizeuy = (int)(maxu/DELTA_U)+1;
    int sizeux = 2*(int)(maxu/DELTA_U);
    
	power = create2DMatrixFloat(nbins,Nfreq);
    residcube = create3DMatrixFloat(sizeuy,sizeux,Nfreq);
    powercube = create3DMatrixFloat(sizeuy,sizeux,Nfreq);
	powertemp = create2DMatrixFloat(nbins,Nfreq);
	powerflag = create2DMatrixFloat(nbins,Nfreq);
	power_resid = create2DMatrixFloat(nbins,Nfreq);
	power_residimag = create2DMatrixFloat(nbins,Nfreq);
	powertot = create2DMatrixFloat(nbins,Nfreq);
	fg_num = create2DMatrixInt(nbins,Nfreq);
	lssa_weights = create2DMatrixFloat(nbins,Nfreq);
    lssa_weights_cube = create3DMatrixFloat(sizeuy,sizeux,Nfreq);
    
    

Tsys = 110.;

    //float Tsys = 440.;
	float kboltz = 1380.;
	float Aeff = 21.;
	float chanwidth = global_chanwidth;
	float time_int = tsample;
    
    if (band == 1)    obsvolume = 1./((2.*DELTAUU*0.363497)*(2.*DELTAUU*0.363497))*chanwidth*Nchan;    // sr Hz   high-band
    if (band == 0)    obsvolume = 1./((2.*DELTAUU*(152./182.)*0.363497)*(2.*DELTAUU*(152./182.)*0.363497))*chanwidth*Nchan;    // sr Hz   low-band
    if ((band != 0)&&(band != 1)) obsvolume = 1./((2.*DELTAUU*(82./182.)*0.363497)*(2.*DELTAUU*(82./182.)*0.363497))*chanwidth*Nchan;    // sr Hz - ultralow band


	float sigma = Tsys*kboltz/Aeff/sqrt(2.*chanwidth*time_int);
//	float sigma=1.;
	sigma2 = sigma*sigma;


/**************** Start OpenMP to parallelize ********************/

omp_set_num_threads(nthreads);
#pragma omp parallel shared (powertot,power,powerflag,power_resid,power_residimag,lssa_weights,space,fg_num,lssa_weights_cube,powercube,residcube)
{

#pragma omp for schedule(dynamic) nowait

/************************************************/
	/* Loop over kbins */
 
  	for (k=0;k<nbins-4;k++){
//    for (k=0;k<20;k++){
        printf("bin number: %d\n",k);
	fprintf(flog,"Starting bin number: %d\n",k);

	int i_in,j_in,k_in,kk,ii,dist_bin;
        unsigned long counter=0;
		float pwr=0.;
	float *totpower=NULL,*crosspower=NULL,*flagpower=NULL,*residpower=NULL,*crossvar=NULL,*residpowerimag=NULL,*norm2=NULL;

	totpower = calloc(Nfreq,sizeof(float));
	crosspower = calloc(Nfreq,sizeof(float));
	flagpower = calloc(Nfreq,sizeof(float));
	residpower = calloc(Nfreq,sizeof(float));
	residpowerimag = calloc(Nfreq,sizeof(float));
	crossvar = calloc(Nfreq,sizeof(float));
	norm2 = calloc(nbins+1, sizeof(float));



	for (i_in=0;i_in<u_size+1;i_in++){
  //      printf("%d\n",i_in);

		for (j_in=0;j_in<2.*u_size;j_in++){

		// determine which bin u,v cell falls within
		for (ii=0;ii<nbins+1;ii++) norm2[ii] = (((uu[i_in]*uu[i_in]+vv[j_in]*vv[j_in]) - lperp[ii]*lperp[ii])*((uu[i_in]*uu[i_in]+vv[j_in]*vv[j_in]) - lperp[ii]*lperp[ii]));

		dist_bin = (int) (sqrt(uu[i_in]*uu[i_in]+vv[j_in]*vv[j_in])/space);
		k_in = minimum_location(norm2,nbins);
			
				//			pwr = 0.;
				//			for (ii=0;ii<Nchan;ii++) pwr += cabs(vis1[i_in][j_in][ii])/((float)Nchan);
			
				//			if (pwr > 20.) printf("Power exceeds max\n");
			
				//		if (pwr < 20.){
			if ((k_in == k)){
 //               printf("Correct bin\n");

			counter = 0;
			for (ii=0;ii<Nchan;ii++){
				if (weights0[i_in][j_in][ii] > 0.) counter++;
			}

            if (counter > Nchan) fprintf(flog,"counter %d\n",counter);
                
			if (counter > Nchan/2){

			  if (dft_cross_pro(lperp[k]+maxu/((float) nbins)/2.,Nchan,vis0[i_in][j_in],vis1[i_in][j_in],flags0[i_in][j_in],flags1[i_in][j_in],weights0[i_in][j_in],totpower,crosspower,flagpower,residpower,residpowerimag,crossvar,sigma2,bias_mode,band,tsample) == 0){

float sum=0.;       
				for (kk=0;kk<Nfreq;kk++){

					if (crossvar[kk] != 0.){
                        

					powertot[k][kk] += (totpower[kk]*crossvar[kk]);
					power[k][kk] += (crosspower[kk]*crossvar[kk]);
                	    powerflag[k][kk] += (flagpower[kk]*crossvar[kk]);
					power_resid[k][kk] += (residpower[kk]*crossvar[kk]);
					power_residimag[k][kk] += (residpowerimag[kk]*crossvar[kk]);
                	    powercube[i_in][j_in][kk] +=(crosspower[kk]*crossvar[kk]);
                	    residcube[i_in][j_in][kk] +=(residpower[kk]*crossvar[kk]);
					lssa_weights_cube[i_in][j_in][kk] += (crossvar[kk]);
		                   lssa_weights[k][kk] += crossvar[kk];
/*					printf("Var: %d %g\n",kk,crossvar[kk]);
					printf("Tot: %d %g\n",kk,totpower[kk]);
					sum += cabs(totpower[kk]);
					printf("sum: %d %g\n",kk,sum); */
					fg_num[k][kk]++;
								
		
					}
		
		
				}

				memset(totpower,0,Nfreq*sizeof(float));
				memset(crossvar,0,Nfreq*sizeof(float));
				memset(flagpower,0,Nfreq*sizeof(float));
				memset(residpower,0,Nfreq*sizeof(float));
				memset(residpowerimag,0,Nfreq*sizeof(float));
				memset(crosspower,0,Nfreq*sizeof(float));
				memset(norm2,0,(nbins+1)*sizeof(float));


			}
			}
					//		}
				
			}

		}  /* end loop over rows */
	
	fflush(flog);
     
	}  //  ******** end loop over columns ********

	free(totpower);
	free(crosspower);
	free(residpower);
	free(crossvar);
	free(flagpower);
	free(residpowerimag);
        free(norm2);

 //       fprintf(flog,"eta bin 0 %f eta bin 1 %f eta bin 10 %f eta bin 20 %f eta bin 30 %f eta bin 100 %f\n",power[k][0],power[k][1],power[k][10],power[k][20],power[k][30],power[k][100]);
	fprintf(flog,"Completed bin number: %d\n",k);
        
	}

}  /*************** END OPENMP PARALLEL LOOP ********************/

    /* write-out files */
	printf("Writing files\n");
    
    sprintf(outfilename1,"%scrosspower_%s_%d.iter.%s.dat",getenv("OUTPUTDIR"),pol,bias_mode,date);
    fprintf(flog,"Cross power file: %s\n",outfilename1);
printf("Cross power file: %s\n",outfilename1);
    if ((fptr = fopen(outfilename1,"w")) == NULL) printf("failed to open file for writing\n");
 //   for (i=0;i<601;i++){
  //      for (j=0;j<1200;j++) fwrite(power[i][j],sizeof(float),Nfreq,fptr);
  //  }
    for (i=0;i<nbins;i++) fwrite(power[i],sizeof(float),Nfreq,fptr);
    fclose(fptr);


    sprintf(outfilename2,"%stotpower_%s_%d.iter.%s.dat",getenv("OUTPUTDIR"),pol,bias_mode,date);
    fptr = fopen(outfilename2,"w");
    for (i=0;i<nbins;i++) fwrite(powertot[i],sizeof(float),Nfreq,fptr);
    fclose(fptr);
    
    sprintf(outfilename3,"%sflagpower_%s_%d.iter.%s.dat",getenv("OUTPUTDIR"),pol,bias_mode,date);
    fptr = fopen(outfilename3,"w");
    for (i=0;i<nbins;i++) fwrite(powerflag[i],sizeof(float),Nfreq,fptr);
    fclose(fptr);
    
    sprintf(outfilename4,"%sresidpower_%s_%d.iter.%s.dat",getenv("OUTPUTDIR"),pol,bias_mode,date);
    fptr = fopen(outfilename4,"w");
    for (i=0;i<nbins;i++) fwrite(power_resid[i],sizeof(float),Nfreq,fptr);
    fclose(fptr);

    sprintf(outfilename6,"%sresidpowerimag_%s_%d.iter.%s.dat",getenv("OUTPUTDIR"),pol,bias_mode,date);
    fptr = fopen(outfilename6,"w");
    for (i=0;i<nbins;i++) fwrite(power_residimag[i],sizeof(float),Nfreq,fptr);
    fclose(fptr);    

    sprintf(outfilename5,"%soutputweights_%s_%d.iter.%s.dat",getenv("OUTPUTDIR"),pol,bias_mode,date);
    fptr = fopen(outfilename5,"w");
  //  for (i=0;i<601;i++){
  //      for (j=0;j<1200;j++) fwrite(lssa_weights[i][j],sizeof(float),Nfreq,fptr);
  //  }
    for (i=0;i<nbins;i++) fwrite(lssa_weights[i],sizeof(float),Nfreq,fptr);
    fclose(fptr);

    sprintf(outfilename5_1,"%sfg_num_%s_%d.iter.%s.dat",getenv("OUTPUTDIR"),pol,bias_mode,date);
    fptr = fopen(outfilename5_1,"w");
     for (i=0;i<nbins;i++) fwrite(fg_num[i],sizeof(int),Nfreq,fptr);
    fclose(fptr);

    
	sprintf(lperpfilename,"%slperp.dat",getenv("OUTPUTDIR"));
	fptr = fopen(lperpfilename,"w");
	fwrite(lperp,sizeof(float),nbins+1,fptr);
	fclose(fptr);

    sprintf(paramsfilename,"%sparams.%s.txt",getenv("OUTPUTDIR"),date);
    fptr = fopen(paramsfilename,"w");
    fprintf(fptr,"%s %f","Volume (sr Hz): ",obsvolume);
    fclose(fptr);
    
	 
    fprintf(flog,"Output files written\n");
    
    fclose(flog);

		/*************** FREE ARRAYS *********************/



return 0;

}




/******************************************************************************************************/


// Functions

/*************************/
// Name: dft_cross_pro
// Function: Computes the LSSA of the spectral information given the sampling in the data
// Date: April 29, 2019
// Author: Trott
/*************************/

int dft_cross_pro(float u, int Nchan, float complex* vis0, float complex* vis1, float complex* flags0, float complex* flags1, float* weights0, float* totpower, float* crosspower, float* flagpower, float* residpower, float* residpowerimag, float* crossvar, float sigma2, int bias_mode, int band, float t_sample){

	int limit=Nchan/2,i,j,counter=0,num=0,M,Nfreq,one,loc=0,freqend,endval,Nchanend,counter2=0,loc2=0;
    FILE *fptr=NULL;
	time_t current_time;
	double *nut_func=NULL;
    double complex alpha,zero,alphaneg,vis0tot,vis1tot,flags0tot,flags1tot;
	int *freq=NULL,*freqref=NULL, *eta=NULL, scalar=128, lwork, iter;
	double complex *C0dbl=NULL , *work=NULL, *work2=NULL, *H=NULL, *Hdag=NULL;
	double complex *HCS0=NULL, *HCS1=NULL;
	double complex *cov0=NULL, *HCF0=NULL, *HCF1=NULL, *F1=NULL, *F0=NULL, *S1=NULL, *S0=NULL;
	float *HCS0_interp=NULL,*HCS1_interp=NULL,*HCF0_interp=NULL,*HCF1_interp=NULL,*HCS2_interp=NULL,*HCS3_interp=NULL;
        float *HCS0_temp=NULL,*HCS1_temp=NULL,*HCF0_temp=NULL,*HCF1_temp=NULL,*HCS2_temp=NULL,*HCS3_temp=NULL;
    float complex *HCS0_shift=NULL,*HCS1_shift=NULL,*HCF0_shift=NULL,*HCF1_shift=NULL,*C0=NULL;
	float *var0=NULL, *var0_shift=NULL, *etaout=NULL, *var0_interp=NULL, *eta_actual=NULL,*normal=NULL, *input_freq=NULL, kernel=0.,freq_term=0,filter_total=0.,bp=1.;
		double a_nut[4]={0.3635819,0.4891775,0.1365995,0.0106411};
	double a_nut2[7]={0.27122036,-0.4334446,0.21800412,-0.06578534,0.010761867,-7.700127e-4,1.36808831e-5};
	double a_nut5[5]={0.3232153788877,0.4714921439576,0.175534129960197,0.00284969901061499,0.0001261357088293};
	double a_nut3[9]={0.2384331152777,-0.400554539,0.235824253,-0.095279188584,0.025373955167,-4.152432907e-3,3.68560416e-4,-1.384355594e-5,1.1618083589e-7};
	char *trans="N",*filetemp[1024],*filetemp2[1024];
	float nut_norm = 0.,t_gs=253.,eta_gs=0.01,u0_gs=10.;

   
	
    a_nut[0]=0.36358192677076108;
	a_nut[1]=-0.489177437145017;
	a_nut[2]=0.136599513978692;
	a_nut[3]=-0.01064112210553003;   // B-N
	
    
    //    a_nut[4]={0.35875,0.48829,0.14128,0.01168};   // B-H

/*
    a_nut[0]=0.35875;
    a_nut[1]=-0.48829;
    a_nut[2]=0.14128;
    a_nut[3]=-0.01168;   // B-H
*/


     
//        a_nut[0]=0.42659;
 //       a_nut[1]=0.49656;
//        a_nut[2]=0.076849;
//        a_nut[3]=0.0;   // Blackman
    
    
	/*
     a_nut[0] = 25./46.;
     a_nut[1] = -21./46.;
     a_nut[2] = 0.;
     a_nut[3] = 0.;
	*/
    // Hamming
    
    
	/*
    a_nut[0] = 1.;
       a_nut[1] = 0.;
      a_nut[2] = 0.;
        a_nut[3] = 0.;  // Rectangle
	*/

	alpha = 1.0;
    alphaneg = -1.0;
	zero = 0.0;
	one = 1;
	input_freq = calloc(Nchan,sizeof(float));
	if (band == 1){
		for (i=0;i<Nchan;i++) input_freq[i] = LOWER_FREQ_HIGH + (float) i*80000.;
	} else {
		for (i=0;i<Nchan;i++) input_freq[i] = LOWER_FREQ + (float) i*80000.;
	}

	nut_func = calloc(Nchan,sizeof(double));
//	for (j=0;j<Nchan;j++){
//		for (i=0;i<4;i++) nut_func[j] += a_nut[i]*cos(2.*M_PI*(float) i*(j-Nchan/2.)/(float) Nchan);
//	}
	for (j=0;j<Nchan;j++){
			//		for (i=0;i<4;i++) nut_func[j] += a_nut[i]*cos(2.*M_PI*(float) i*((float) j)/((float) Nchan-1.));
		for (i=0;i<4;i++) nut_func[j] += a_nut[i]*cos(2.*M_PI*(float) i*((float) j)/((float) Nchan-1.));
		//		printf("nut func %f\n",nut_func[j]);
	}
	

	freqref = calloc(Nchan,sizeof(int));
	for (i=0;i<Nchan;i++) freqref[i] = i;

//	printf("Made it here 1\n");

    /***************************************************************************************************************/
	/* Set-up spectral coverage and LSSA size for this u,v */
    /***************************************************************************************************************/

	filter_total=0.;
	float in_sq1 = 0.;
	float in_sq0 = 0.;

	for (i=0;i<Nchan;i++){
		in_sq1 += (vis1[i])*conj(vis1[i]);
        in_sq0 += (vis0[i])*conj(vis0[i]);
		if (weights0[i] != 0.){
			 num++;
			filter_total += sqrt(nut_func[i]);
		}
	}
    
//printf("Taper normalisation %f\n",filter_total);

	// return failure if there are too few spectral points for LSSA
    if (num < limit){
        
        printf("Too few spectral points %d\n",num);
    return 1;
    }
    

	Nfreq = Nchan;
    num = (int)((float)Nchan*0.875);
//    num = Nchan;
 //   if (bias_mode == 10) num=166;
    if ((bias_mode == 15)||(bias_mode == 14)) num = 162;
	M = num;
    
 //   printf("Channel numbers: %d %d %d %d\n",Nchan,Nfreq,M,num);

//	if (M != Nchan) printf("Incorrect sizes Nchan %d M %d num %d\n",Nchan,M,num);
    
    /***************************************************************************************************************/
	/* Create full covariance matrix in freq-freq space using sum of thermal weights and FG term */
    /***************************************************************************************************************/


    /* Add foreground covariance to radiometric covariance, if necessary */
    
    C0 = calloc(Nchan,sizeof(float complex));
    
    counter=0;
    counter2=0;

    float weighttot = 0;
    for (i=0;i<Nchan;i++) weighttot += weights0[i]/(float)Nchan;
    
  //  printf("sigma %f weighttot %f\n",sigma2,weighttot);
    
    float weight2tot = 0.;
    for (i=0;i<Nchan;i++) weight2tot += sqrt(weights0[i])/(float)num;
    
  /*  for (i=0;i<Nchan;i++){
        
        if (weights0[i] > 0.) weighttot += 1./((float)Nchan*weights0[i]);
    }
    weighttot = 1./weighttot;
*/
 
//    printf("Weights going in: %f %f\n",weighttot,weight2tot);
    
	for (i=0;i<Nchan;i++){
  //      printf("weights %f\n",weights0[i]);
		if (weights0[i] != 0){
                    loc = (int) (counter2);

//	weights0[i] = 1.;

	C0[loc] = sigma2/weighttot;

            counter2++;
            
             //          printf("counter2 %d C0 %f\n",counter2,creal(C0[loc]));
    }
	}

    
	//	printf("Made it here\n");


    /***************************************************************************************************************/
	/* Factorize full and FG covariance matrices, and then invert for LSSA */
    /***************************************************************************************************************/

    

    /***************************************************************************************************************/
	/* Define visibility and flags vectors and operate on them with C^-1 */
    /***************************************************************************************************************/


	S0 = calloc(num,sizeof(double complex));
	F0 = calloc(num,sizeof(double complex));
	S1 = calloc(num,sizeof(double complex));
	F1 = calloc(num,sizeof(double complex));
	freq = calloc(num,sizeof(int));
  
  
    vis0tot = 0.;
    vis1tot = 0.;
    flags0tot = 0.;
    flags1tot = 0.;
    
    
if (flag_dc){

printf("Flagging DC mode\n");
//fprintf(flog,"Flagging DC mode\n");

    for (i=0;i<M;i++){
        vis0tot += vis0[i]/M;
        vis1tot += vis1[i]/M;
        flags0tot += flags0[i]/M;
        flags1tot += flags1[i]/M;
    }
    
}


	counter=0;
	for (i=0;i<Nchan;i++){
		if (weights0[i] != 0){
            		freq[counter] = freqref[i];

            S0[counter] = (vis0[i]-vis0tot)*nut_func[i]/(filter_total)*(float)num*sqrt(weights0[i])/weight2tot; // This is  S = (BdagB)^-1 Bdag V
            F0[counter] = (flags0[i]-flags0tot)*nut_func[i]/filter_total*(float)num*sqrt(weights0[i])/weight2tot;
            S1[counter] = (vis1[i]-vis1tot)*nut_func[i]/(filter_total)*(float)num*sqrt(weights0[i])/weight2tot; // This is S
            F1[counter] = (flags1[i]-flags1tot)*nut_func[i]/filter_total*(float)num*sqrt(weights0[i])/weight2tot;
 
			counter++;
		}
	}

    
    /***************************************************************************************************************/
	/* Define Fourier basis functions for LSSA - H and H^dagger */
    /***************************************************************************************************************/

  
	// shift frequency
	int tempfreq = freq[0];
	for (i=0;i<M;i++){
	 freq[i] -= tempfreq;
	}
    

	eta = malloc(Nfreq*sizeof(int));
	eta_actual = malloc(Nfreq*sizeof(float));

	for (i=0;i<Nfreq;i++){
		eta[i] = i;
//		eta_actual[i] = (float)i/(float)Nfreq*(float)Nchan/2.;
        if (Nfreq % 2 == 0){ eta_actual[i] = ((float)i/(float)Nfreq - 0.5 + 1./(float)Nfreq)*(float)Nchan; } else { eta_actual[i] = ((float)i/(float)Nfreq - 0.5 + 1./2./(float)Nfreq)*(float)Nchan; }
//        printf("i %d eta %f\n",i,eta_actual[i]);
	}
 
    
      freqend = freq[num-1];  // last sampled frequency
      Nchanend = freqend+1;
 
    Nchanend = Nchan;
    
    if (Nchanend != Nchan) printf("Nchanend %d\n",Nchanend);
    

	// define spectral matrix
	H = malloc(M*Nfreq*sizeof(double complex));
	Hdag = malloc(M*Nfreq*sizeof(double complex));
 
	// compute H matrices and noise matrix

	for (j=0;j<M;j++){
		for (i=0;i<Nfreq;i++){
			H[j*Nfreq + i] = cexp(I*(2.*M_PI*(double)freq[j]*(double)eta[i]/(double)Nfreq))/sqrt((float)M);
		}
        
	}

        
    for (i=0;i<Nfreq;i++){
		for (j=0;j<M;j++) Hdag[j + i*M] = cexp(I*(-2.*M_PI*(double)freq[j]*(double)eta[i]/(double)Nfreq))/sqrt((float)M);
	}
    

    /***************************************************************************************************************/
	/* Perform matrix multiplication for full covariance matrix: (H^T C^-1 H) */
    /***************************************************************************************************************/
	
	cov0 = calloc(Nfreq,sizeof(double complex));

	weighttot = 0;     
    		for (i=0;i<Nchan;i++) weighttot += C0[i];
		for (i=0;i<Nfreq;i++) cov0[i] = ((float) Nfreq)/weighttot;

 //   printf("weighttot %f\n",weighttot);

    /***************************************************************************************************************/
	/* Compute spectral modes of data: Hdagger C^-1 S */
    /***************************************************************************************************************/

	HCS0 = calloc(Nfreq,sizeof(double complex));
	cblas_zgemv(CblasRowMajor,CblasNoTrans,Nfreq,M,&alpha,Hdag,M,S0,1,&zero,HCS0,1 );

	HCF0 = calloc(Nfreq,sizeof(double complex));
	cblas_zgemv(CblasRowMajor,CblasNoTrans,Nfreq,M,&alpha,Hdag,M,F0,1,&zero,HCF0,1 );

	HCS1 = calloc(Nfreq,sizeof(double complex));
	cblas_zgemv(CblasRowMajor,CblasNoTrans,Nfreq,M,&alpha,Hdag,M,S1,1,&zero,HCS1,1 );

	HCF1 = calloc(Nfreq,sizeof(double complex));
	cblas_zgemv(CblasRowMajor,CblasNoTrans,Nfreq,M,&alpha,Hdag,M,F1,1,&zero,HCF1,1 );


    /***************************************************************************************************************/
	/* Compute the product of the covariance matrices: Cinv^dagger # Cinv = cov0^T # cov0 */
    /***************************************************************************************************************/

    //   printf("Computing error bars...\n");
    var0 = calloc(Nfreq,sizeof(float));
    for (i=0;i<Nfreq;i++) var0[i] = creal(cov0[i]*cov0[i])/16.;
    
     
    /***********************************************************************************************************/
 
    /***************************************************************************************************************/
	/* Shift modes to run from -ve to +ve spectral modes, and interpolate output LSSA and weights to regular kpar grid */
    /***************************************************************************************************************/

   
	/* Interpolate the outputs for each onto the regular eta grid */

	etaout = calloc(Nchan,sizeof(float));
	for (i=0;i<Nchan;i++) etaout[i] = (int) ((float)i-Nchan/2 + 1.);

	HCS0_interp = calloc(Nfreq,sizeof(float));
	HCS1_interp = calloc(Nfreq,sizeof(float));
	HCF0_interp = calloc(Nfreq,sizeof(float));
	HCF1_interp = calloc(Nfreq,sizeof(float));
	HCS2_interp = calloc(Nfreq,sizeof(float));
	HCS3_interp = calloc(Nfreq,sizeof(float));

	var0_interp = calloc(Nchan,sizeof(float));

	HCS0_temp = calloc(Nchan,sizeof(float));
	HCS1_temp = calloc(Nchan,sizeof(float));
	HCF0_temp = calloc(Nchan,sizeof(float));
	HCF1_temp = calloc(Nchan,sizeof(float));
	HCS2_temp = calloc(Nchan,sizeof(float));
	HCS3_temp = calloc(Nchan,sizeof(float));    

    // Shift frequencies to be sequential
    
	HCS0_shift = calloc(Nfreq,sizeof(float complex));
	HCS1_shift = calloc(Nfreq,sizeof(float complex));
	HCF0_shift = calloc(Nfreq,sizeof(float complex));
	HCF1_shift = calloc(Nfreq,sizeof(float complex));
	var0_shift = calloc(Nfreq,sizeof(float));
    
    if (Nfreq % 2 == 0){
        endval = (int) Nfreq/2.;
    //    printf("Even\n");
        
        for (i=0;i<Nfreq/2.-1;i++){
            HCS0_shift[i] = HCS0[(int) (i + Nfreq/2.) + 1];
            HCS1_shift[i] = HCS1[(int) (i + Nfreq/2.) + 1];
            HCF0_shift[i] = HCF0[(int) (i + Nfreq/2.) + 1];
            HCF1_shift[i] = HCF1[(int) (i + Nfreq/2.) + 1];
            var0_shift[i] = var0[(int) (i + Nfreq/2.) + 1];
        }
        for (i=Nfreq/2.-1;i<Nfreq;i++){
            HCS0_shift[i] = HCS0[(int) (i - Nfreq/2.) + 1];
            HCS1_shift[i] = HCS1[(int) (i - Nfreq/2.) + 1];
            HCF0_shift[i] = HCF0[(int) (i - Nfreq/2.) + 1];
            HCF1_shift[i] = HCF1[(int) (i - Nfreq/2.) + 1];
            var0_shift[i] = var0[(int) (i - Nfreq/2.) + 1];
         }
        
    } else {
        endval = (int) Nfreq/2. + 1;
   //     printf("Odd\n");

        for (i=0;i<(int) Nfreq/2.;i++){
            HCS0_shift[i] = HCS0[i + (int) (Nfreq/2.)+1];
            HCS1_shift[i] = HCS1[i + (int) (Nfreq/2.)+1];
            HCF0_shift[i] = HCF0[i + (int) (Nfreq/2.)+1];
            HCF1_shift[i] = HCF1[i + (int) (Nfreq/2.)+1];
            var0_shift[i] = var0[i + (int) (Nfreq/2.)+1];
        }
        for (i=(int) Nfreq/2.;i<Nfreq;i++){
            HCS0_shift[i] = HCS0[i - (int) (Nfreq/2.)];
            HCS1_shift[i] = HCS1[i - (int) (Nfreq/2.)];
            HCF0_shift[i] = HCF0[i - (int) (Nfreq/2.)];
            HCF1_shift[i] = HCF1[i - (int) (Nfreq/2.)];
            var0_shift[i] = var0[i - (int) (Nfreq/2.)];
        }

    
    }
 
 
 /*   current_time = time(NULL);
    sprintf(filetemp,"/data/gridded_vis/covcov_%f.dat",u);
    fptr = fopen(filetemp,"w");
    fwrite(cov0cov0,sizeof(double complex),Nfreq*Nfreq,fptr);
    fclose(fptr);
 
    sprintf(filetemp,"/data/gridded_vis/even%d_%f.dat",current_time,u);
    fptr = fopen(filetemp,"w");
    fwrite(HCS1_shift,sizeof(float complex),Nfreq,fptr);
    fclose(fptr);
   */
    
 //   printf("Interpolating...\n");

	/* form power and interpolate */

	for (i=0;i<Nfreq;i++){
		HCS0_interp[i] = HCS0_shift[i]*conj(HCS0_shift[i]);
		HCF0_interp[i] = HCF0_shift[i]*conj(HCF0_shift[i]);
		HCS1_interp[i] = HCS1_shift[i]*conj(HCS1_shift[i]);
		HCF1_interp[i] = HCF1_shift[i]*conj(HCF1_shift[i]);
		HCS2_interp[i] = creal(HCS0_shift[i]*conj(HCS1_shift[i]) + HCS1_shift[i]*conj(HCS0_shift[i]));
		HCS3_interp[i] = cimag(HCS0_shift[i]*conj(HCS1_shift[i]) + HCS1_shift[i]*conj(HCS0_shift[i]));
	}


    /***************************************************************************************************************/
	/* Compute outputs */
    /***************************************************************************************************************/


	for (i=0;i<Nfreq;i++){
		totpower[i] = (HCS1_interp[i]);

        residpower[i] = 0.25*((HCS1_shift[i]+HCS0_shift[i])*conj(HCS1_shift[i]+HCS0_shift[i]) - (HCS1_shift[i]-HCS0_shift[i])*conj(HCS1_shift[i]-HCS0_shift[i]));

    }
    

    for (i=0;i<Nfreq;i++){
		residpowerimag[i] = (HCS3_interp[i]);
        crosspower[i] = 0.25*(totpower[i] - (HCS0_interp[i]));
		flagpower[i] = 0.25*(HCF1_interp[i] - HCF0_interp[i]);
		crossvar[i] = fabs(var0_shift[i]);
 //       printf("%f\n",crossvar[i]);
	}
    
//	printf("In power tot: %f	In power res: %f	Out power tot: %f	Out power resid: %f	Ratio: %f\n",in_sq1,in_sq0,out_sq1,out_sq0,out_sq1/in_sq1);

//	clean_up();
   
	if (etaout != NULL) free(etaout);
	if (eta_actual != NULL) free(eta_actual);
	if (eta != NULL) free(eta);
	if (freq != NULL) free(freq);
	if (freqref != NULL) free(freqref);
	if (cov0 != NULL) free(cov0);
    if (HCS0 != NULL) free(HCS0);
	if (HCS1 != NULL) free(HCS1);
	if (HCF0 != NULL) free(HCF0);
	if (HCF1 != NULL) free(HCF1);
	if (HCS0_interp != NULL) free(HCS0_interp);
	if (HCS2_interp != NULL) free(HCS2_interp);
	if (HCS3_interp != NULL) free(HCS3_interp);
	if (HCS1_interp != NULL) free(HCS1_interp);
	if (HCF0_interp != NULL) free(HCF0_interp);
	if (HCF1_interp != NULL) free(HCF1_interp);
	if (HCS0_temp != NULL) free(HCS0_temp);
	if (HCS1_temp != NULL) free(HCS1_temp);
	if (HCS2_temp != NULL) free(HCS2_temp);
	if (HCS3_temp != NULL) free(HCS3_temp);
	if (HCF0_temp != NULL) free(HCF0_temp);
	if (HCF1_temp != NULL) free(HCF1_temp);
    if (HCS0_shift != NULL) free(HCS0_shift);
	if (HCS1_shift != NULL) free(HCS1_shift);
	if (HCF0_shift != NULL) free(HCF0_shift);
	if (HCF1_shift != NULL) free(HCF1_shift);
	if (var0 != NULL) free(var0);
    if (var0_shift != NULL) free(var0_shift);
	if (var0_interp != NULL) free(var0_interp);
//	if (R != NULL) free(R);

	if (H != NULL) free(H);
	if (Hdag != NULL) free(Hdag);
    if (S0 != NULL) free(S0);
    if (S1 != NULL) free(S1);
    if (F0 != NULL) free(F0);
    if (F1 != NULL) free(F1);
//	if (cov_fg != NULL) free_memFloat(cov_fg,Nchan,Nchan);

   
	return 0;

}






void interlin(float complex* input, int sizein, float* input_loc, float* output_loc, int sizeout, float complex* output) {
	/* interpolate input array at output_loc locations and place into output */
	int i,loclow;
	float loc_des,delta,loc,frac;
	float complex m=0,c=0;
	float sumin=0.,sumout=0;;

//	for (i=0;i<sizein;i++) sumin += cabs(input[i]);

	delta = input_loc[1] - input_loc[0];
	//c = input[0];
  //  printf("sizein %d sizeout %d delta %f inputloc0 %f c %f\n",sizein,sizeout,delta,input_loc[0],creal(c));
    
	// loop through desired locations
	for (i=0;i<sizeout;i++){
        
		loc_des = output_loc[i];
        
		if (loc_des <= input_loc[0]){
			output[i] = (input[0]);
		} else if (loc_des >= input_loc[sizein-1]){
			output[i] = (input[sizein-1]);
		} else {
			loc = (loc_des - input_loc[0])/delta;  // location in input array of desired point
			loclow = (int)loc;
			frac = ((float) (loc_des - input_loc[loclow]))/delta;
			m = (input[loclow+1]-input[loclow]);
			c = input[loclow];
			output[i] = m*(frac) + c;
	printf("i %d loc %f loc_des %f input_loc[loclow] %f input[loclow] %f input[loclow+1] %f output[i] %f power %f\n",i,loc,loc_des,input_loc[loclow],cimag(input[loclow]),cimag(input[loclow+1]),cimag(output[i]),output[i]*conj(output[i]));
  //          printf("i %d m %f c %f loclow %d frac %f loc %f\n",i,creal(m),creal(c),loclow,frac,loc);
		}
 //       printf("i %d %f %f %f %f %f\n",i,loc_des,creal(output[i]),cimag(output[i]),creal(m),creal(c));

	}
 //   printf("Completed\n");

//	for (i=0;i<sizeout;i++) sumout += cabs(output[i]);

//	printf("Sumin: %f\n",sumin);
//	printf("Sumout: %f\n",sumout);
//
//	for (i=0;i<sizeout;i++) output[i] = output[i]*sumin/sumout;

}


void interlinf(float* input, int sizein, float* input_loc, float* output_loc, int sizeout, float* output) {
	/* interpolate input array at output_loc locations and place into output */
	int i,loclow;
	float loc_des,delta,loc,frac;
	float m,c;
	float sumin=0.,sumout=0;;

	for (i=0;i<sizein;i++) sumin += (input[i]);

	delta = input_loc[1] - input_loc[0];
	//c = input[0];
  //  printf("sizein %d sizeout %d delta %f inputloc0 %f c %f\n",sizein,sizeout,delta,input_loc[0],creal(c));
    
	// loop through desired locations
	for (i=0;i<sizeout;i++){
        
		loc_des = output_loc[i];
        
		if (loc_des <= input_loc[0]){
			output[i] = (input[0]);
		} else if (loc_des >= input_loc[sizein-1]){
			output[i] = (input[sizein-1]);
		} else {
			loc = (loc_des - input_loc[0])/delta;  // location in input array of desired point
			loclow = (int)loc;
			frac = ((float) (loc_des - input_loc[loclow]))/delta;
			m = (input[loclow+1]-input[loclow]);
			c = input[loclow];
			output[i] = m*(frac) + c;
  //          printf("i %d loc %f loc_des %f input_loc[loclow] %f input[loclow] %f input[loclow+1] %f output[i] %f\n",i,loc,loc_des,input_loc[loclow],input[loclow],input[loclow+1],output[i]);
		}
 //       printf("i %d %f %f %f %f %f\n",i,loc_des,creal(output[i]),cimag(output[i]),creal(m),creal(c));

	}
 //   printf("Completed\n");

	for (i=0;i<sizeout;i++) sumout += (output[i]);

//	printf("Sumin: %f\n",sumin);
//	printf("Sumout: %f\n",sumout);

	for (i=0;i<sizeout;i++) output[i] = output[i]*sumin/sumout;

}


void interpad(float complex* input, int sizein, float* input_loc, float* output_loc, int sizeout, float complex* output) {
	/* interpolate input array at output_loc locations and place into output */
	int i,j;
	float loc_des,tol=1.e-2;
    
	// loop through desired locations
	for (i=0;i<sizeout;i++){
        
		loc_des = output_loc[i];
        
		if (loc_des <= input_loc[0]){
			output[i] = (input[0]);
		} else if (loc_des >= input_loc[sizein-1]){
			output[i] = (input[sizein-1]);
		} else {
		for (j=0;j<sizein;j++){			
			if ((loc_des-input_loc[j])*(loc_des-input_loc[j]) < tol){ 
				output[i] = input[j];
		//		printf("Match!\n");
		}
	//		printf("j %d i %d input_loc[j] %f loc_des %f\n",j,i,input_loc[j],loc_des);
		}
            
		}
 //       printf("i %d %f %f %f %f %f\n",i,loc_des,creal(output[i]),cimag(output[i]),creal(m),creal(c));

	}
 //   printf("Completed\n");

}


void interpadf(float* input, int sizein, float* input_loc, float* output_loc, int sizeout, float* output) {
	/* interpolate input array at output_loc locations and place into output */
	int i,loclow,j;
	float loc_des,tol=1.e-2;

	// loop through desired locations
	for (i=0;i<sizeout;i++){
        
		loc_des = output_loc[i];
        
		if (loc_des <= input_loc[0]){
			output[i] = (input[0]);
		} else if (loc_des >= input_loc[sizein-1]){
			output[i] = (input[sizein-1]);
		} else {
		for (j=0;j<sizein;j++){			
			if ((loc_des-input_loc[j])*(loc_des-input_loc[j]) < tol) output[i] = input[j];
		}
  //          printf("i %d m %f c %f loclow %d frac %f loc %f\n",i,creal(m),creal(c),loclow,frac,loc);
		}
 //       printf("i %d %f %f %f %f %f\n",i,loc_des,creal(output[i]),cimag(output[i]),creal(m),creal(c));

	}
 //   printf("Completed\n");

}



double bessj1( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          1 at input x                                      */
/*------------------------------------------------------------*/
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
         +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
      ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
         +y*(99447.43394+y*(376.9991397+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-2.356194491;
      ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
         +y*(0.2457520174e-5+y*(-0.240337019e-6))));
      ans2=0.04687499995+y*(-0.2002690873e-3
         +y*(0.8449199096e-5+y*(-0.88228987e-6
         +y*0.105787412e-6)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
      if (x < 0.0) ans = -ans;
   }
   return ans;
}


int minimum_location(float* array, int size) {
    int c, location = 1;
	float *minimum;
 
    minimum = array;
    *minimum = *array;
 
    for ( c = 1 ; c < size ; c++ ) 
    {
        if ( *(array+c) < *minimum ) 
        {
           *minimum = *(array+c);
           location = c+1;
        }
    } 
 
    return location;
}



double trace_mat(double **matrix, int size){
  int i;
  double out=0.;

  for (i=0;i<size;i++) out += matrix[i][i];

  //  printf("out %g\n",out);
  return out;


}

/*
int round(float num){

int n = (int)(num < 0 ? (num - 0.5) : (num + 0.5));

 return n;
}
*/

void where(long *vec_in, long int size, int limit, long int *vec_out, int *counter){
  int i,c=0;

  for (i=0;i<size;i++){
    if (vec_in[i] == limit){
      vec_out[c] = i;
      c++;
    }
  }

  counter = &c;
  //     printf("c %i\n",c);

}

int matrixVectorMultiply(int vsize, double vec[], double **mat, double outVector[]){
int i,j;

for (i=0;i<vsize;i++){
	outVector[i] = 0.;
	for (j=0;j<vsize;j++){
		outVector[i] += vec[j]*mat[i][j];
		}
}
return 0;
}

int matrixMatrixMultiply(int x1size, int y1size, int x2size, int y2size, double **mat1, double **mat2, double **outMat){
int i,j,k;

if (x1size != y2size) exit(1);

for (i=0;i<x2size;i++){
for (j=0;j<y1size;j++){
	for (k=0;k<x1size;k++){

		outMat[i][j] += mat1[k][j]*mat2[i][k];
	
	}
}
}
return 0;
}

int CmatrixMatrixMultiply(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag){
int i,j,k;

 if (x1size != y2size) exit(1);
 
 for (i=0;i<x2size;i++){
   for (j=0;j<y1size;j++){
     for (k=0;k<x1size;k++){
	  
       outMatreal[i][j] += mat1real[k][j]*mat2real[i][k] - mat1imag[k][j]*mat2imag[i][k];
       outMatimag[i][j] += mat1real[k][j]*mat2imag[i][k] + mat1imag[k][j]*mat2real[i][k];
	  
     }
   }
 }
 return 0;
 
}

int CmatrixMatrixMultiplyConj(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag){
int i,j,k;

/* Conjugates the imaginary part of the second matrix */

if (x1size != y2size) exit(1);

for (i=0;i<x2size;i++){
for (j=0;j<y1size;j++){
	for (k=0;k<x1size;k++){

		outMatreal[i][j] += mat1real[k][j]*mat2real[i][k] + mat1imag[k][j]*mat2imag[i][k];
		outMatimag[i][j] += -mat1real[k][j]*mat2imag[i][k] + mat1imag[k][j]*mat2real[i][k];
	
	}
}
}
return 0;

}

int CmatrixMatrixMultiplyConjTrans(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag){
int i,j,k;

/* Conjugates the imaginary part of the second matrix, and transposes */

/* Must be symmetric!! */

if (x1size != y2size) exit(1);

for (i=0;i<x2size;i++){
for (j=0;j<y1size;j++){
	for (k=0;k<x1size;k++){

		outMatreal[i][j] += mat1real[k][j]*mat2real[k][i] + mat1imag[k][j]*mat2imag[k][i];
		outMatimag[i][j] += -mat1real[k][j]*mat2imag[k][i] + mat1imag[k][j]*mat2real[k][i];
	
	}
}
}
return 0;

}

int matrixTranspose(int xsize, int ysize, double **mat_in, double **mat_out){
int i,j;

for (i=0;i<xsize;i++){
for (j=0;j<ysize;j++){

	mat_out[j][i] = mat_in[i][j];

}
}
return 0;
}

double vectorVectorMultiply(int vsize, double vec1[], double vec2[]){
int i;
double temp;

	temp = 0.;
	for (i=0;i<vsize;i++){
		temp += vec1[i]*vec2[i];
	}
	return temp;

}

double **create2DMatrix(int sizex, int sizey){
double **output;
int k;

	output = calloc(sizex,sizeof(double*));
	for (k=0;k<sizex;k++){
		output[k] = calloc(sizey,sizeof(double));
	}
	return output;

}


float complex **create2DMatrixComplex(int sizex, int sizey){
float complex **output;
int k;

	output = calloc(sizex,sizeof(float complex*));
	for (k=0;k<sizex;k++){
		output[k] = calloc(sizey,sizeof(float complex));
	}
	return output;

}



double ***create3DMatrix(int sizex, int sizey, int sizez){
double ***output;
int k,l;

	output = calloc(sizex,sizeof(double**));
	for (k=0;k<sizex;k++){
		output[k] = calloc(sizey,sizeof(double*));
	}
	for (k=0;k<sizex;k++){
	for (l=0;l<sizey;l++){
		output[k][l] = calloc(sizez,sizeof(double));
	}
	}
	return output;

}

/*
float ***create3DMatrixFloat(int sizex, int sizey, int sizez){
float ***output;
int k,l;

	output = calloc(sizex,sizeof(float**));
	for (k=0;k<sizex;k++){
		output[k] = calloc(sizey,sizeof(float*));
	}
	for (k=0;k<sizex;k++){
	for (l=0;l<sizey;l++){
		output[k][l] = calloc(sizez,sizeof(float));
	}
	}
	return output;

}
*/

float complex ***create3DMatrixComplex(int sizex, int sizey, int sizez){
float complex ***output;
int k,l;

	output = calloc(sizex,sizeof(float complex**));
	for (k=0;k<sizex;k++){
		output[k] = calloc(sizey,sizeof(float complex*));
	}
	for (k=0;k<sizex;k++){
	for (l=0;l<sizey;l++){
		output[k][l] = calloc(sizez,sizeof(float complex));
	}
	}
	return output;

}



int ***create3DMatrixInt(int sizex, int sizey, int sizez){
int ***output;
int k,l;

	output = calloc(sizex,sizeof(int**));
	for (k=0;k<sizex;k++){
		output[k] = calloc(sizey,sizeof(int*));
	}
	for (k=0;k<sizex;k++){
	for (l=0;l<sizey;l++){
		output[k][l] = calloc(sizez,sizeof(int));
	}
	}
	return output;

}

long ***create3DMatrixLong(int sizex, int sizey, int sizez){
long ***output;
int k,l;

	output = calloc(sizex,sizeof(long**));
	for (k=0;k<sizex;k++){
		output[k] = calloc(sizey,sizeof(long*));
	}
	for (k=0;k<sizex;k++){
	for (l=0;l<sizey;l++){
		output[k][l] = calloc(sizez,sizeof(long));
	}
	}
	return output;

}

double ****create4DMatrix(int sizex, int sizey, int sizez, int sizew){
double ****output;
 int k,l,m;

 output = calloc(sizex,sizeof(double***));
 for (k=0;k<sizex;k++){
   output[k] = calloc(sizey,sizeof(double**));
 }
 for (k=0;k<sizex;k++){
   for (l=0;l<sizey;l++){
     output[k][l] = calloc(sizez,sizeof(double*));
   }
 }
 for (k=0;k<sizex;k++){
   for (l=0;l<sizey;l++){
     for (m=0;m<sizez;m++){
       output[k][l][m] = calloc(sizew,sizeof(double));
     }
   }
 }
 return output;


}




long ****create4DMatrixLong(int sizex, int sizey, int sizez, int sizew){
long ****output;
 int k,l,m;

 output = calloc(sizex,sizeof(long***));
 for (k=0;k<sizex;k++){
   output[k] = calloc(sizey,sizeof(long**));
 }
 for (k=0;k<sizex;k++){
   for (l=0;l<sizey;l++){
     output[k][l] = calloc(sizez,sizeof(long*));
   }
 }
 for (k=0;k<sizex;k++){
   for (l=0;l<sizey;l++){
     for (m=0;m<sizez;m++){
       output[k][l][m] = calloc(sizew,sizeof(long));
     }
   }
 }
 return output;


}


float **create2DMatrixFloat(int sizex, int sizey){
float **output;
int k;

	output = calloc(sizex,sizeof(float*));
	for (k=0;k<sizex;k++){
		output[k] = calloc(sizey,sizeof(float));
	}
	return output;

}

int **create2DMatrixInt(int sizex, int sizey){
int **output;
int k;

	output = calloc(sizex,sizeof(int*));
	for (k=0;k<sizex;k++){
		output[k] = calloc(sizey,sizeof(int));
	}
	return output;

}

long int **create2DMatrixLong(int sizex, int sizey){
long int **output;
int k;

	output = calloc(sizex,sizeof(long int*));
	for (k=0;k<sizex;k++){
		output[k] = calloc(sizey,sizeof(long int));
	}
	return output;

}

unsigned long **create2DMatrixShort(int sizex, int sizey){
unsigned long **output;
int k;

	output = calloc(sizex,sizeof(unsigned long*));
	for (k=0;k<sizex;k++){
		output[k] = calloc(sizey,sizeof(unsigned long));
	}
	return output;

}

double *create1DVector(int size){
double *output;
	
	output = calloc(size,sizeof(double));
	return output;

}


void memzero2D(double **matrix,int xsize, int ysize){
  int i;

  for (i=0;i<xsize;i++) memset(matrix[i],0.,ysize*sizeof(double));

}

void memzero2DFloat(float **matrix,int xsize, int ysize){
    int i;
    
    for (i=0;i<xsize;i++) memset(matrix[i],0.,ysize*sizeof(float));
    
}


void memzero2Dint(int **matrix,int xsize, int ysize){
  int i;

  for (i=0;i<xsize;i++) memset(matrix[i],0.,ysize*sizeof(int));

}

void memzero2Dlong(long int **matrix,int xsize, int ysize){
  int i;

  for (i=0;i<xsize;i++) memset(matrix[i],0.,ysize*sizeof(long int));

}


void memzero3D(double ***matrix,int xsize, int ysize, int zsize){
  int i,j;

  for (i=0;i<xsize;i++){
    for (j=0;j<ysize;j++) memset(matrix[i][j],0.,zsize*sizeof(double));
  }

}


void memzero3Dint(int ***matrix,int xsize, int ysize, int zsize){
  int i,j;

  for (i=0;i<xsize;i++){
    for (j=0;j<ysize;j++) memset(matrix[i][j],0.,zsize*sizeof(int));
  }

}

void free_mem(double **matrix, int xsize, int ysize){
  int i;

  for (i=0;i<xsize;i++) free(matrix[i]);

}

void free_memFloat(float **matrix, int xsize, int ysize){
    int i;
    
    for (i=0;i<xsize;i++) free(matrix[i]);
    
}


void free_mem3D(double ***matrix, int xsize, int ysize, int zsize){
    int i,j;
    
    for (i=0;i<xsize;i++){
        for (j=0;j<ysize;j++) free(matrix[i][j]);
            }

}




