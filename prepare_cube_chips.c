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
#include "uvfits.h"
#include "primary_beam.h"
#include "fitsio.h"
#include <omp.h>

// gcc -Wall -g prepare_lssa.c cspline.c uvfits.c -o prepare_lssa -L/usr/local/lib -lm -lcfitsio

#define expand_factor 1

/* Function prototypes */
void free_mem3D(double ***matrix, int xsize, int ysize, int zsize);
int omp_get_num_threads(void);
double vectorVectorMultiply(int vsize, double vec1[], double vec2[]);
void free_memFloat(float **matrix, int xsize, int ysize);
int matrixVectorMultiply(int vsize, double vec[], double **mat, double outVector[]);
double **create2DMatrix(int sizex, int sizey);
double ***create3DMatrix(int sizex, int sizey, int sizez);
int ***create3DMatrixInt(int sizex, int sizey, int sizez);
long ***create3DMatrixLong(int sizex, int sizey, int sizez);
double ****create4DMatrix(int sizex, int sizey, int sizez, int sizew);
void memzero2DFloat(float **matrix,int xsize, int ysize);
long ****create4DMatrixLong(int sizex, int sizey, int sizez, int sizew);
double complex **create2DMatrixComplex(int sizex, int sizey);
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

/* global variables */
float global_period = PERIOD;
float global_umax = UMAX;
float global_chanwidth = CHAN_WIDTH;
float global_lower_freq = LOWER_FREQ_HIGH;

int debug=0;
char *infilename=NULL;
FILE *fpd;

void usage() {
    fprintf(stderr,"Usage:\n");
    fprintf(stderr,"test_readuvfits <options> -i inputfilename\n");
    fprintf(stderr,"\t-d debug_level (0=no debugging)\n");
    exit(0);
}

void parse_cmdline(const int argc,char * const argv[]) {
    int result=0;
    const char optstring[] = "i:p:u:n:c:d:";

    while ( (result = getopt(argc, argv, optstring)) != -1 ) {
        switch (result) {
          case 'i': infilename = optarg;
            break;
            case 'p': global_period = atof(optarg);
                break;
			case 'u': global_umax = atof(optarg);
				break;
            case 'n': global_lower_freq = atof(optarg);
                break;
            case 'c': global_chanwidth = atof(optarg);
                printf("Setting global chanwidth to: %f\n",global_chanwidth);
                break;
          case 'd': debug = atoi(optarg);
            break;
          default:
              fprintf(stderr,"unknown option: %c\n",result);
              usage(argv[0]);
        }
    }
}





/* Main */
int main (int argc, char *argv[]) {

  int i,j,ch, size_fin, freq_index_start,set, freq_index, Nchan, size_half, band;
  FILE *fptrv=NULL,*fpbb=NULL,*flog=NULL;
  float frequency=LOWER_FREQ,*freq=NULL, *weights=NULL,*norm=NULL,*w1=NULL,LOWER_FRE=0.,zero_offset = 0.,weight_tot=0.,weight_max=0.,umx=300.;
  float tot_tot = 0.,diff_tot=0.,diff_max=0.,tot_max=0.,mx_w=0.;
  double chanwidth=0.,freq_start;
  double u_size;
  char *date=NULL,syslogfile[1024],*pol=NULL,bvfilename[1024],bbfilename[1024],normfilename[1024],filename_real[1024],filename_real2[1024],filename_real6[1024],*ext=NULL;
  unsigned long counter=0;
  float complex *diff=NULL,*tot=NULL, *outfold2=NULL;
  double complex **outw=NULL,*outfold=NULL,**out=NULL,**outdiff=NULL;
    float spec=1.6,factor_scale=1.;
	
  fpd = stderr;

// Information


if (argc < 2){
    fprintf(stderr,"Usage: <options> %s input_ext Nchan freq_index_start pol ext band\n",argv[0]);
    fprintf(stderr,"\t -p period (seconds)\n");
    fprintf(stderr,"\t -c chanwidth\n");
	fprintf(stderr,"\t -u umax\n");
    fprintf(stderr,"\t -n bottom frequency (Hz)\n");
    exit(1);
}


printf("Program to combine data over frequency.\n");

date = argv[1];
Nchan = atoi(argv[2]);
freq_index_start = atoi(argv[3]);

pol = argv[4];
ext = argv[5];
band = atoi(argv[6]);

    parse_cmdline(argc,argv);
    
    printf("Global chanwidth: %f chanwidth %f\n",global_chanwidth,chanwidth);
    
    chanwidth = global_chanwidth;
    LOWER_FRE = global_lower_freq;
	
	umx = global_umax;
    
    printf("LOWER_FREQ_HIGH %f umax %f\n",LOWER_FRE,umx);

    
    printf("Band: %d Chanwidth %f Period %f\n",band,chanwidth,global_period);

/* Open log file */
	sprintf(syslogfile,"%ssyslog_preparediff_%s_%d.txt",getenv("OUTPUTDIR"),date,freq_index_start);
	printf("syslogfile: %s\n",syslogfile);
	if ((flog=fopen(syslogfile,"w")) == NULL){
	fprintf(stderr,"Cannot open output log file\n");
	return 1;
	}

    time_t tim = time(NULL);
    struct tm *tmptr = gmtime(&tim);

fprintf(flog,"Processed with version %f of CHeIIPS DIFF on %s\n",VERSION,asctime(localtime(&tim)));

printf("Global chanwidth: %f chanwidth %f global_period %f Global umax: %f\n",global_chanwidth,chanwidth,global_period,global_umax);

    /* define uv grid size and spacing */

u_size = floor(umx/DELTA_U);

if (HALF_PLANE_FLAG){size_fin = (int) 2.*u_size*u_size;} else { size_fin = (int) 4.*u_size*u_size;}

size_half = (int) 2.*u_size*(u_size+1);
    
    printf("u size %g size_half %d size_fin %d\n",u_size,size_half,size_fin);

freq = calloc(Nchan,sizeof(float));
if (band == 1){
	for (i=0;i<Nchan;i++) freq[i] = (LOWER_FRE + freq_index_start*COARSE_CHAN_WIDTH + chanwidth*i)/1.e6;
} else {
	for (i=0;i<Nchan;i++) freq[i] = (LOWER_FREQ + freq_index_start*COARSE_CHAN_WIDTH + chanwidth*i)/1.e6;
}

    printf("LOWER_FREQ_HIGH %f\n",LOWER_FRE);
    
    /**************************************************/
    
    weights = calloc(size_half*Nchan,sizeof(float));
    w1 = calloc(size_half*Nchan,sizeof(float));
    outw = create2DMatrixComplex(size_fin,Nchan);
    outfold = calloc(size_half*Nchan,sizeof(float complex));
    

    /* Loop over channels to compute each frequency */
    /* Compute weights contribution */
    
    sprintf(bvfilename,"%sweightsc_%s.%s.dat",getenv("OUTPUTDIR"),pol,date);
    
    printf("File: %s\n",bvfilename);
    
    /* Compute bv contribution */
    
    if ((fptrv = fopen(bvfilename,"r")) == NULL){
        // file does not exist - ignore
        printf("File does not exist %s... ignoring...\n",bvfilename);
        
    } else {
        
        printf("filename: %s %d %d\n",bvfilename,size_fin,Nchan);
        
        mx_w = 0.;
        
        for (i=0;i<size_fin;i++){
            for (j=0;j<Nchan;j++){
                    // if (fread((outw[i]),sizeof(double complex),Nchan,fptrv) != Nchan) fprintf(flog,"Error: input file is the incorrect size. File: %s\n",bvfilename);
                fread(&outw[i][j],sizeof(double complex),1,fptrv);
                if (creal(outw[i][j]) > mx_w){
                    mx_w=creal(outw[i][j]);
                }
            }
        }
        
        printf("Weights max: %f\n",mx_w);
        
        fclose(fptrv);
        
    }
    
    printf("%f %f\n",creal(outw[320400][400]),creal(outw[200000][300]));
    
	for (ch=0;ch<Nchan;ch++){
        //	for (ch=3;ch<4;ch++){  //lower 8 MHz
        // get frequency information //
        
        frequency = freq[ch];

   //     factor_scale = pow((freq[Nchan-1])/(freq[ch]),spec);
    //    printf("frequency %f chan %d factor_scale %f ratio %f\n",frequency,ch,factor_scale,(freq[Nchan-1])/(freq[ch]));
        
	//	fprintf(flog,"pol %s\n",pol);
        
   //     fprintf(flog,"frequency %g\n",frequency);
   //     printf("frequency %g\n",frequency);
        
            
   //         printf("Folding\n");
            
            /* Fold onto half-plane */
            for (i=1;i<u_size+1;i++){
                for (j=0;j<2*u_size;j++){
                    
                    int index1 = i*(2.*u_size) + j;
                    int index2 = (float)(u_size*2-1-i)*(float)(u_size*2) + u_size*2-1-j;
                    int index3 = (float)(2*u_size*(u_size+1)) -1 - index1;
              //      printf("%d %d %d\n",index1,index2,index3);
                //    printf("%d %d %d %f %f\n",index1,index2,index3,creal(outw[index1][ch]),creal(outw[index2][ch]));
                    outfold[index3] = ((outw[index1][ch]) + conj(outw[index2][ch]));
                }
            }
            
            
     //       printf("Moving\n");
            
            
            // annoying row-major order
            for (i=0;i<size_half;i++){
                if (creal(outfold[i]) > 0.){
            //
                    w1[i*Nchan + ch] = creal(outfold[i]);///norm[ch];
               //     printf("i %d weight %f\n",i,w1[i*Nchan+ch]);
                }
            }
            
        
        
          /* Form final weights */
        
        for (i=0;i<size_half;i++){
            if ((w1[i*Nchan+ch] != 0.)){
                weights[i*Nchan + ch] = (w1[i*Nchan+ch]);
             } else {
                    weights[i*Nchan + ch] = 0.;
                  }
                
            }
            
            
            
        }
        
        // write_out data for visibility contributions
        
        sprintf(filename_real2,"%sweights_%s.%s.dat",getenv("OUTPUTDIR"),pol,ext);
        
        fptrv = fopen(filename_real2,"w");
        if(fptrv==NULL) {
            printf("Error: can't open output file of vis.\n");
        } else {
            printf("File opened successfully.\n");
            for (i=0;i<size_half*Nchan;i++){
                //      double temp=creal(out_array[i]);
		if (weights[i] > weight_max) weight_max = weights[i];
		weight_tot += weights[i];
                fwrite(&((weights[i])),sizeof(float),1,fptrv);
            }
        }
        fclose(fptrv);
    
    

/* Set-up output arrays */

diff = calloc(size_half*Nchan,sizeof(float complex));
tot = calloc(size_half*Nchan,sizeof(float complex));
out = create2DMatrixComplex(size_fin,Nchan);
    outdiff = create2DMatrixComplex(size_fin,Nchan);
    outfold2 = calloc(size_half*Nchan,sizeof(float complex));
//outfold = calloc(size_half*Nchan,sizeof(float complex));


    
/************************************************/
/* Loop over channels to compute each frequency */
 
    sprintf(bvfilename,"%sbv_%s.%s.dat",getenv("OUTPUTDIR"),pol,date);
    sprintf(bbfilename,"%sbvdiff_%s.%s.dat",getenv("OUTPUTDIR"),pol,date);

    printf("File: %s\n",bvfilename);

    /* Compute bv contribution */

    if ((fptrv = fopen(bvfilename,"r")) == NULL){
            // file does not exist - ignore
        fprintf(flog,"File does not exist %s... ignoring...\n",bvfilename);
        
    }

fprintf(flog,"filename: %s\n",bvfilename);

    for (i=0;i<size_fin;i++){
        for (j=0;j<Nchan;j++){
            fread(&out[i][j],sizeof(double complex),1,fptrv);
        }
    }
    
fclose(fptrv);
    
    
    
        if ((fptrv = fopen(bbfilename,"r")) == NULL){
            // file does not exist - ignore
            fprintf(flog,"File does not exist %s... ignoring...\n",bbfilename);

        } else {
            
            fprintf(flog,"filename: %s\n",bbfilename);
            
            for (i=0;i<size_fin;i++){
                for (j=0;j<Nchan;j++){
                    fread(&outdiff[i][j],sizeof(double complex),1,fptrv);
                }
            }
            
            fclose(fptrv);
            
        }

	for (ch=0;ch<Nchan;ch++){
//	for (ch=3;ch<4;ch++){  //lower 8 MHz
	// get frequency information //

	frequency = freq[ch];


	//	fprintf(flog,"pol %s\n",pol);

  //  fprintf(flog,"frequency %g\n",frequency);
	//		printf("frequency %g\n",frequency);


       //     printf("Folding\n");
            
            /* Fold onto half-plane */
            for (i=1;i<u_size+1;i++){
                for (j=0;j<2*u_size;j++){
                    
                    int index1 = i*(2.*u_size) + j;
                    int index2 = (float)(u_size*2-1-i)*(float)(u_size*2) + u_size*2-1-j;
                    int index3 = (float)(2*u_size*(u_size+1)) -1 - index1;
                    outfold[index3] = (out[index1][ch] + conj(out[index2][ch]));
                        //           if (ch > 0)   printf("outfold %lf out out %lf %lf\n",creal(outfold2[index3]),creal(out[index1]),creal(out[index2]));
                }
            }
            
            
      //      printf("Moving\n");
            
            
                // annoying row-major order
            for (i=0;i<size_half;i++){
                    //		printf("%f %f %f\n",creal(outfold[i]),w1[i*Nchan+ch],w2[i*Nchan+ch]);
                if ((w1[i*Nchan + ch] != 0.)){
                        //			if (ch > 0) printf("weights %f %d %d\n",w1[i*Nchan+ch],i*Nchan+ch,i);
                    tot[i*Nchan + ch] = outfold[i]/(w1[i*Nchan+ch]*factor_scale + zero_offset);
                        //		if (ch > 0) printf("totals %f %f %f %d\n",creal(outfold2[i]),cimag(tot[i*Nchan + ch]),w1[i*Nchan+ch],ch);
                }
            }
            
        

				// read-in second set data
			//	printf("Second set\n");

	

			//		printf("Folding\n");

					/* Fold onto half-plane */
					for (i=1;i<u_size+1;i++){
						for (j=0;j<2*u_size;j++){

							int index1 = i*(2.*u_size) + j;
							int index2 = (float)(u_size*2-1-i)*(float)(u_size*2) + u_size*2-1-j;
							int index3 = (float)(2*u_size*(u_size+1)) -1 - index1;
							outfold2[index3] = (outdiff[index1][ch] + conj(outdiff[index2][ch]));
						}
					}

			//		printf("Moving\n");


					// annoying row-major order
					for (i=0;i<size_half;i++){
					  if ((w1[i*Nchan + ch] != 0.)){
                            diff[i*Nchan + ch] = outfold2[i]/(w1[i*Nchan+ch]*factor_scale + zero_offset);
                        }
					}


        

		}
        

		// write_out data for visibility contributions

		sprintf(filename_real2,"%svis_diff_%s.%s.dat",getenv("OUTPUTDIR"),pol,ext);
		sprintf(filename_real,"%svis_tot_%s.%s.dat",getenv("OUTPUTDIR"),pol,ext);

			fptrv = fopen(filename_real2,"w");
			if(fptrv==NULL) {
				printf("Error: can't open output file of vis.\n");
			} else {
			printf("File opened successfully.\n");
			for (i=0;i<size_half*Nchan;i++){
     			//      double temp=creal(out_array[i]);

			if (creal(diff[i]) > diff_max) diff_max = creal(diff[i]);
			diff_tot += creal(diff[i]);
				fwrite(&(diff[i]),sizeof(float complex),1,fptrv);
   				}
			}
			fclose(fptrv);

			fptrv = fopen(filename_real,"w");
			if(fptrv==NULL) {
				printf("Error: can't open output file of vis.\n");
			} else {
			printf("File opened successfully.\n");
			for (i=0;i<size_half*Nchan;i++){
     			//      double temp=creal(out_array[i]);
			if (creal(tot[i]) > tot_max) tot_max = creal(tot[i]);
			tot_tot += creal(tot[i]);
				fwrite(&(tot[i]),sizeof(float complex),1,fptrv);
   				}
			}
			fclose(fptrv);


		/* write completion to log */
		fprintf(flog,"Completed frequency %g, pol %s\n",frequency,pol);

	
	fflush(flog);

    
    /************************************************/
    /* Loop over channels to compute each frequency */
    /* Compute noise contribution */
    
    memset(diff,0.,size_half*Nchan*sizeof(float complex));
    memset(tot,0.,size_half*Nchan*sizeof(float complex));
    
    sprintf(bvfilename,"%snoisec_%s.%s.dat",getenv("OUTPUTDIR"),pol,date);
    sprintf(bbfilename,"%snoisecdiff_%s.%s.dat",getenv("OUTPUTDIR"),pol,date);
    
    printf("File: %s\n",bvfilename);
    
    /* Compute bv contribution */

    if ((fptrv = fopen(bvfilename,"r")) == NULL){
            // file does not exist - ignore
        fprintf(flog,"File does not exist %s... ignoring...\n",bvfilename);
        
    }

fprintf(flog,"filename: %s\n",bvfilename);

    for (i=0;i<size_fin;i++){
        for (j=0;j<Nchan;j++){
            fread(&out[i][j],sizeof(double complex),1,fptrv);
        }
    }
    
fclose(fptrv);
    
    
    
        if ((fptrv = fopen(bbfilename,"r")) == NULL){
            // file does not exist - ignore
            fprintf(flog,"File does not exist %s... ignoring...\n",bbfilename);

        } else {
            
            fprintf(flog,"filename: %s\n",bbfilename);
            
            for (i=0;i<size_fin;i++){
                for (j=0;j<Nchan;j++){
                    fread(&outdiff[i][j],sizeof(double complex),1,fptrv);
                }
            }
            
            fclose(fptrv);
            
            
        }
    
	for (ch=0;ch<Nchan;ch++){
        //	for (ch=3;ch<4;ch++){  //lower 8 MHz
        // get frequency information //
        
        frequency = freq[ch];
        
		fprintf(flog,"pol %s\n",pol);
        
        fprintf(flog,"frequency %g\n",frequency);
        printf("frequency %g\n",frequency);
        
            
    //        printf("Folding\n");
            
            /* Fold onto half-plane */
            for (i=1;i<u_size+1;i++){
                for (j=0;j<2*u_size;j++){
                    
                    int index1 = i*(2.*u_size) + j;
                    int index2 = (float)(u_size*2-1-i)*(float)(u_size*2) + u_size*2-1-j;
                    int index3 = (float)(2*u_size*(u_size+1)) -1 - index1;
                    outfold2[index3] = out[index1][ch] + conj(out[index2][ch]);
                }
            }
            
            
       //     printf("Moving\n");
            
            
            // annoying row-major order
            for (i=0;i<size_half;i++){
	      if ((w1[i*Nchan + ch] != 0.)){
                    tot[i*Nchan + ch] = outfold2[i]/(w1[i*Nchan+ch]*factor_scale+zero_offset);
                 }
            }
            
        
        
        // read-in second set data
 

            
         //   printf("Folding\n");
            
            /* Fold onto half-plane */
            for (i=1;i<u_size+1;i++){
                for (j=0;j<2*u_size;j++){
                    
                    int index1 = i*(2.*u_size) + j;
                    int index2 = (float)(u_size*2-1-i)*(float)(u_size*2) + u_size*2-1-j;
                    int index3 = (float)(2*u_size*(u_size+1)) -1 - index1;
                    outfold2[index3] = outdiff[index1][ch] + conj(out[index2][ch]);
                }
            }
            
       //     printf("Moving\n");
            
            
            // annoying row-major order
            for (i=0;i<size_half;i++){
	      if ((w1[i*Nchan + ch] != 0.)){
                    diff[i*Nchan + ch] = outfold2[i]/(w1[i*Nchan+ch]*factor_scale+zero_offset);
                }
            }
            
            
        
        
    }
    
    
    // write_out data for visibility contributions
    
    sprintf(filename_real2,"%snoise_diff_%s.%s.dat",getenv("OUTPUTDIR"),pol,ext);
    sprintf(filename_real,"%snoise_tot_%s.%s.dat",getenv("OUTPUTDIR"),pol,ext);
    
    fptrv = fopen(filename_real2,"w");
    if(fptrv==NULL) {
        printf("Error: can't open output file of vis.\n");
    } else {
        printf("File opened successfully.\n");
        for (i=0;i<size_half*Nchan;i++){
            //      double temp=creal(out_array[i]);
            fwrite(&(diff[i]),sizeof(float complex),1,fptrv);
        }
    }
    fclose(fptrv);
    
    fptrv = fopen(filename_real,"w");
    if(fptrv==NULL) {
        printf("Error: can't open output file of vis.\n");
    } else {
        printf("File opened successfully.\n");
        for (i=0;i<size_half*Nchan;i++){
            //      double temp=creal(out_array[i]);
            fwrite(&(tot[i]),sizeof(float complex),1,fptrv);
        }
    }
    fclose(fptrv);
    
    
    /* write completion to log */
    fprintf(flog,"Completed frequency %g, pol %s\n",frequency,pol);

    
    
    fclose(flog);

		/*************** FREE ARRAYS *********************/



return 0;

}





// Functions

double complex **create2DMatrixComplex(int sizex, int sizey){
double complex **output;
int k;

    output = calloc(sizex,sizeof(double complex*));
    for (k=0;k<sizex;k++){
        output[k] = calloc(sizey,sizeof(double complex));
    }
    return output;

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




