#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
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

// gcc -Wall -g grid_vis.c cspline.c uvfits.c -o gridvis -L/usr/local/lib -lm -lcfitsio -lcholmod


/* Function prototypes */
void free_mem3D(double ***matrix, int xsize, int ysize, int zsize);
int omp_get_num_threads(void);
double vectorVectorMultiply(int vsize, double vec1[], double vec2[]);
void free_memFloat(float **matrix, int xsize, int ysize);
void free_memComplex(double complex **matrix, int xsize, int ysize);
int matrixVectorMultiply(int vsize, double vec[], double **mat, double outVector[]);
double gauss_noise (double value1);
void memzero3Dc(double complex ***matrix,int xsize, int ysize, int zsize);
double **create2DMatrix(int sizex, int sizey);
double ***create3DMatrix(int sizex, int sizey, int sizez);
int ***create3DMatrixInt(int sizex, int sizey, int sizez);
long ***create3DMatrixLong(int sizex, int sizey, int sizez);
double ****create4DMatrix(int sizex, int sizey, int sizez, int sizew);
void memzero2DFloat(float **matrix,int xsize, int ysize);
long ****create4DMatrixLong(int sizex, int sizey, int sizez, int sizew);
//float **create2DMatrixFloat(int sizex, int sizey);
unsigned long **create2DMatrixShort(int sizex, int sizey);
int **create2DMatrixInt(int sizex, int sizey);
long int **create2DMatrixLong(int sizex, int sizey);
int matrixMatrixMultiply(int x1size, int y1size, int x2size, int y2size, double **mat1, double **mat2, double **outMat);
int CmatrixMatrixMultiply(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag);
int CmatrixMatrixMultiplyConj(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag);
int CmatrixMatrixMultiplyConjTrans(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag);
int matrixTranspose(int xsize, int ysize, double **mat_in, double **mat_out);
double complex **create2DMatrixComplex(int sizex, int sizey);
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
float global_chanwidth = CHAN_WIDTH;
float global_lower_freq = LOWER_FREQ;
float global_umax = UMAX;
int field = FIELD;
int debug=0;
char *infilename=NULL;
int ewflag=0;
FILE *fpd;

void usage() {
    fprintf(stderr,"Usage:\n");
    fprintf(stderr,"test_readuvfits <options> -i inputfilename\n");
    fprintf(stderr,"\t-d debug_level (0=no debugging)\n");
    exit(0);
}


void parse_cmdline(const int argc,char * const argv[]) {
    int result=0;
    const char optstring[] = "i:e:p:c:f:n:d:u:";

    while ( (result = getopt(argc, argv, optstring)) != -1 ) {
        switch (result) {
          case 'i': infilename = optarg;
            break;
        case 'e': ewflag = 1;
                break;
          case 'p': global_period = atof(optarg);
            break;
          case 'c': global_chanwidth = atof(optarg);
            break;
			case 'u': global_umax = atof(optarg);
			break;
        case 'f': field = atoi(optarg);
            break;
            case 'n': global_lower_freq = atof(optarg);
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

  int size_fin, flag_tsys=0,bandchan,band,bandlowchan,i,j,ch,secondset=0;
  FILE *fpobs=NULL,*ftsysobs=NULL,*flog=NULL,*fraobs=NULL,*flut=NULL;
  float lmst,ha_rad,dec_rad,alt_point,az_point,ha_point,l_rot,m_rot,umx=300.,vis_ratio=1.;
	int pol,polindex,polarray[2]={0,1};
  double ra_app,dec_app,ra_point,dec_point,ra_offset,dec_offset,ra_field,dec_field,ra_zenith,dec_zenith;
  double tsys[2];
    long total_size=936000000;
  double u_size,mjd=0.0,az,el,systemp;
  char obsidfilename[1024],obsidbandfilename[1024],obsidtsysfilename[1024],syslogfile[1024],*ext=NULL,lutfile14m_ew[1024];
  unsigned long obsid;
  int point=4,nthreads=1;
    int **lut14ew=NULL,*flag14m=NULL;
  beamdata *beam_data;
        int Ntriad=56;
    
    FILE *fptr=NULL,*fptrv=NULL,*fptrflags=NULL;
    int w_index,obs_year,obs_month,obs_day,beam_size,kernel_size;
    double checksum_real=0.0,checksum_imag=0.0,checksum_abs=0.0,checksum_cent=0.0,beam_sigma2=0.05;
    float delt=DELTAUU;
    char filename_real[1024],filename_real2[1024],filename_reald[1024],filename_real2d[1024],filename_flags[1024],debugfile[1024];
    double complex  **bdag_v_real=NULL,**bdag_v_realdiff=NULL,**temp_v_real=NULL,**temp_noise=NULL,**temp_v_realdiff=NULL,**temp_noisediff=NULL,**noise_array=NULL,**noisediff_array=NULL,**flag_array=NULL,**temp_flags=NULL;
        float vis_rtot,vis_rdiff,vis_itot,vis_idiff,w_temp=0.,kernel_size_float,frequency_lower=0.,factor_scale=1.,l8,additup=0.;
    
        double a_nut[4]={0.3635819,0.4891775,0.1365995,0.0106411};
        float nut_norm = 0.;
	
  nthreads = omp_get_max_threads();
  printf("max threads: %d\n",nthreads);

		//	 nthreads=1;


  //  fitsfile *fptr,*fptri;
  int res=0,chunk=0;
  uvdata *data,*data2;
  uvReadContext *iter,*iter2;

  fpd = stderr;

// Information


if (argc < 2){
    fprintf(stderr,"Usage: %s <options> uvfits_filename obs_id extension band \n",argv[0]);
    fprintf(stderr,"\t -p period (seconds)\n");
    fprintf(stderr,"\t -ew flag 14m EW baselines\n");
    fprintf(stderr,"\t -c chanwidth (Hz)\n");
    fprintf(stderr,"\t -n bottom frequency (Hz)\n");
    fprintf(stderr,"\t -field fieldnum (0=EoR0, 1=EoR1, >2 anything)\n");
	fprintf(stderr,"\t -u umax\n");
    exit(1);
}
    
    
printf("Program to calculate the gridded visibilities for an MWA observation.\n\n");

infilename = argv[1];
obsid = atol(argv[2]);
ext = argv[3];
band = atoi(argv[4]);
  //  pol = atoi(argv[5]);

    parse_cmdline(argc,argv);

{

	char tmp1[FILENAME_MAX],*p;


/* Open log file */

    sprintf(syslogfile,"%ssyslog%10ld_all.txt",getenv("OUTPUTDIR"),obsid);
	printf("syslogfile: %s\n",syslogfile);
	if ((flog=fopen(syslogfile,"w")) == NULL){
	fprintf(stderr,"Cannot open output log file\n");
	return 1;
	}
}

    time_t tim = time(NULL);
    struct tm *tmptr = gmtime(&tim);

fprintf(flog,"Processed with version %f of CHIPS on %s\n",VERSION,asctime(localtime(&tim)));

	umx = global_umax;

 
    
    
    
        printf("Global period: %f, Global chanwidth: %f\n",global_period,global_chanwidth);
        

/* define uv grid size and spacing */

u_size = floor(umx/DELTA_U);

if (HALF_PLANE_FLAG){size_fin = (int) 2.*u_size*u_size;} else { size_fin = (int) 4.*u_size*u_size;}

printf("Num. u bins: %lg, size of uv plane: %d\n",u_size,size_fin);

//WMAX = ((NUM_W_PLANES/2)*DELTA_W + (NUM_W_STACK/2)*DELTA_WSTACK)*1.0;
//WMAX = NUM_W_STACK*DELTA_WSTACK;

printf("WMAX %lg\n",WMAX);
        
    
    
    // Loop over polarisations!!!!
    
    for (pol=0;pol<2;pol++){
        
        secondset=0;
        
        /* read-in uvfits data */
        
        uvfitsSetDebugLevel(0);
        res = readUVFITSInitIterator(infilename, &data, &iter);
        
        /* Define dense vector for bdag_v accumulation and flags for uv-sampling */
        
        bdag_v_real = create2DMatrixComplex(size_fin,data->n_freq);
        noise_array = create2DMatrixComplex(size_fin,data->n_freq);
        bdag_v_realdiff = create2DMatrixComplex(size_fin,data->n_freq);
        noisediff_array = create2DMatrixComplex(size_fin,data->n_freq);
        flag_array = create2DMatrixComplex(size_fin,data->n_freq);
        
            // define total size of these
        total_size = (long)((float)size_fin * (float)data->n_freq);
        printf("Total size %ld\n",total_size);
        
        
        /* set-up beam to contain enough entries to comfortably contain beam
         region */
        
        beam_size = (int) BEAM_SIZE_FLOAT;
            //  printf("beam size %d\n",beam_size);
        
        if (((float)beam_size)/2. == (int)((float)(beam_size)/2.)) beam_size = beam_size+1;  /* make odd */
        
        kernel_size_float = beam_size*beam_size;
        kernel_size = (int)kernel_size_float;
        
        
        
        a_nut[0]=0.36358192677076108;
        a_nut[1]=0.489177437145017;
        a_nut[2]=0.136599513978692;
        a_nut[3]=0.01064112210553003;   // B-N
        
        
        
        res = readUVFITSInitIterator(infilename, &data2, &iter2);
        
        if (res !=0) {
            fprintf(stderr,"readUVFITSInitIterator failed on timestep2 with error %d\n",res);
            return res;
        }
        
        fflush(stdout);
        
        /* Dummy read of data2 to increment it to second timestep */
        if ((res=readUVFITSnextIter(data2,iter2)) !=0){
            printf("Second set iteration not found at start\n");
            return 42;
        }
        
        
        /* code runs entirely within WHILE loop, which loops over sections of the input data file (one timestep per iteration) */
        
        uvfitsSetDebugLevel(0);
        printf("freq. channels %d, num vis %d\n",data->n_freq,data->n_vis);
            //printf("uv %d\n",readUVFITSnextIter(data,iter));
        while ((res=readUVFITSnextIter(data,iter)) ==0) {
            fprintf(stdout,"Chunk %d. Time: %f. baselines: %d\n",chunk++,data->date[0],data->n_baselines[0]);
            
            
            if ((res=readUVFITSnextIter(data2,iter2)) !=0){
                printf("Second set iteration not found\n");
                secondset=42;
            }
            
            if (secondset != 42){
                
                
                /************************************************/
                
                /* Collect relevant information from uvfits input file */
                
                printf("Number of frequencies: %d\n",data->n_freq);
                    //   printf("data delta frequency %g data central frequency %g number channels %d\n",data->freq_delta,data->cent_freq,data->n_freq);
                
                /**************** Start OpenMP to parallelize ********************/
                
                omp_set_num_threads(nthreads);
#pragma omp parallel private (ch) shared (a_nut,bdag_v_real, noise_array, bdag_v_realdiff, noisediff_array, flag_array)
                {
                    
                    
                    
                    int *include_vis=NULL,loc=0,flag,flag2,flagbase,xx,top,ii,jj,k,iii,jjj,freq_index;
                    double *beam_real=NULL,*beam_imag=NULL,*u_lex_small=NULL,*v_lex_small=NULL,uu,vv,ww,norm,*beamsq_real=NULL,*beamsq_imag=NULL,*beamsq_real2=NULL,*beamsq_imag2=NULL;
                    double distance=0.,uu2,vv2,ww2,*u_lex_small2=NULL,*v_lex_small2=NULL,*beam_real2=NULL,*beam_imag2=NULL,u_loc1,v_loc1,frequency;
                    float vis_r1,vis_r2,vis_i1,vis_i2,v1,v2,v3,v4,weight,weight_norm,normal,total_norm_sq;
                    double complex inter=0.+I*0.,interorig=0.+I*0.;
                    long contrib_vis=0;
                    
                    
                    
                    
#pragma omp for
                    
                    
                    
                    /************************************************/
                    /* Loop over channels to compute each frequency separately for each OpenMP thread */
                    
                    for (ch=0;ch<data->n_freq;ch++){
                        
                        
                        u_lex_small = calloc(kernel_size,sizeof(double));
                        v_lex_small = calloc(kernel_size,sizeof(double));
                        u_lex_small2 = calloc(kernel_size,sizeof(double));
                        v_lex_small2 = calloc(kernel_size,sizeof(double));
                        
                        beam_real = calloc(kernel_size,sizeof(double));
                        beam_imag = calloc(kernel_size,sizeof(double));
                        
                        beam_real2 = calloc(kernel_size,sizeof(double));
                        beam_imag2 = calloc(kernel_size,sizeof(double));
                        
                        beamsq_real = calloc(kernel_size,sizeof(double));
                        beamsq_imag = calloc(kernel_size,sizeof(double));
                        
                        beamsq_real2 = calloc(kernel_size,sizeof(double));
                        beamsq_imag2 = calloc(kernel_size,sizeof(double));
                        
                        include_vis = calloc(data->n_baselines[0],sizeof(int));
                        
                        
                            // get frequency information //
                        
                        freq_index = (((data->cent_freq)-data->n_freq/2*data->freq_delta)+ch*data->freq_delta-LOWER_FREQ)/COARSE_CHAN_WIDTH;
                        frequency = ((data->cent_freq)-data->n_freq/2*data->freq_delta)+ch*data->freq_delta;
                        
                        
                            //       printf("Frequency: %g baselines: %d\n",frequency,data->n_baselines[0]);
                        /* Loop over the two polarisations */
                        
                            //   for (polindex=0;polindex<NPOL;polindex++){
                            //       for (polindex=0;polindex<1;polindex++){
                        
                            //   pol = polarray[polindex];
                        
                            //	printf("pol %d\n",pol);
                        
                        systemp = 177.;
                        
                        
                        
                        
                        /*********************************************************************************************/
                        /*********************************************************************************************/
                        
                        srand48((long)time(NULL));
                        
                        norm=1.;
                        
                        /* define frequency at low end of coarse channel */
                        
                        frequency_lower = LOWER_FREQ;
                            //  factor_scale = frequency/frequency_lower;
                            //    factor_scale = 1.;
                            //  printf("factor_scale %f\n",factor_scale);
                        
                            //   printf("Here2?\n");
                        /*****************************************************************/
                        /* Loop over baselines - first loop to compute memory allocation */
                        contrib_vis = 0;
                            //    printf("Here?\n");
                            //    for (iii=0;i<data->n_baselines[0];i++){ include_vis[iii]=0; }
                            //     memset(include_vis,0,sizeof(int)*data->n_baselines[0]);
                        
                            //    printf("Here3?\n");
                        
                        
                        for (iii=0;iii<data->n_baselines[0];iii++){
                                //   printf("baseline: %d\n",i);
                            
                            /* Compute whether uv point lies within range of uv-grid and wstack */
                            
                            distance = sqrt(data->u[0][iii]*data->u[0][iii]+data->v[0][iii]*data->v[0][iii])*(frequency);
                                //           printf("i %d distance: %lg w: %g, weight %g, WMAX %f\n",i,distance,data->w[0][iii]*frequency,data->weightdata[0][(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)],WMAX);
                                //  printf("Distance: %f\n",distance);
                            
                            
                            if ((distance >= 0.95*(umx-SIZE_BEAM/2*INTRINSIC_DELTA_U)) || (distance < 0.5) || (sqrt(data->w[0][iii]*frequency*data->w[0][iii]*frequency) >= WMAX)|| (data->weightdata[0][(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] <= 0.)){
                                
                            } else {
                                    //  printf("Included\n");
                                contrib_vis++;
                                include_vis[iii] = 1;
                                
                            }
                            
                        }
                        
                        if (debug) fprintf(flog,"contrib vis: %ld\n",contrib_vis);
                        
                        /*****************************************************************/
                        /* Loop over baselines - second loop to actually process data */
                        
                            //     printf("Contributing vis: %d\n",contrib_vis);
                        
                        if (contrib_vis == 0){
                        } else {
                            
                            
                            for (iii=0;iii<data->n_baselines[0];iii++){
                                    //  printf("baseline: %d\n",i);
                                
                                
                                if (include_vis[iii] != 1){
                                } else {
                                    
                                        //contrib_vis++;
                                    
                                    uu = data->u[0][iii]*(frequency);
                                    vv = data->v[0][iii]*(frequency);
                                    ww = data->w[0][iii]*(frequency);
                                    uu2 = data2->u[0][iii]*(frequency);
                                    vv2 = data2->v[0][iii]*(frequency);
                                    ww2 = data2->w[0][iii]*(frequency);
                                    flag = 0;
                                    
                                        //  printf("basenum %d uu vv ww: %f %f %f\n",iii,uu,vv,ww);
                                    
                                    
                                    if (((HALF_PLANE_FLAG)&&(uu < 0.0))||(ww < 0.)){
                                            // test         if (((HALF_PLANE_FLAG)&&(uu < 0.0))){
                                        flag = 1;
                                        uu = -uu;
                                        vv = -vv;
                                        ww = -ww;
                                        uu2 = -uu2;
                                        vv2 = -vv2;
                                        ww2 = -ww2;
                                        
                                        
                                        
                                    }
                                    
                                    /* remove points in -ve v for u=0 */
                                    if ((HALF_PLANE_FLAG)&&(round(uu/DELTA_U) == 0)&&(vv < 0.)){
                                        flag = 1;
                                        vv = -vv;
                                        vv2 = -vv2;
                                    }
                                    
                                    if (debug) fprintf(flog,"u %g v %g w %g\n",uu,vv,ww);
                                    if(isnan(uu+vv+ww) == 1) printf("u %g v %g w %g\n",uu,vv,ww);
                                    
                                    if (debug) fflush(flog);
                                    
                                    
                                    loc = 0;
                                    
                                    /* lexicographic u and v co-ordinates of small *sky* Fourier plane around the Fourier location */
                                    for (k=0;k<beam_size;k++){
                                        for (jjj=0;jjj<beam_size;jjj++){
                                            
                                            u_lex_small[loc] = (double) ((k)-(int) beam_size/2)*DELTA_U + round(uu/DELTA_U)*DELTA_U;
                                            v_lex_small[loc] = (double) ((jjj)-(int) beam_size/2)*DELTA_U + round(vv/DELTA_U)*DELTA_U;
                                            u_lex_small2[loc] = (double) ((k)-(int) beam_size/2)*DELTA_U + round(uu2/DELTA_U)*DELTA_U;
                                            v_lex_small2[loc] = (double) ((jjj)-(int) beam_size/2)*DELTA_U + round(vv2/DELTA_U)*DELTA_U;
                                            if (debug) printf("loc %ld u_lex %g v_lex %g\n",loc,u_lex_small[loc],v_lex_small[loc]);
                                            loc++;
                                        }
                                    }
                                    
                                    
                                    
                                        //          printf("Here beam\n");
                                    
                                    /* Compute beam gridding by an analytic Blackman-Harris, independent of frequency, w and pointing */
                                    
                                        //    DELTAUU = 5.;
                                    delt = DELTAUU;
                                    interorig = 0. + I*0.;
                                    
                                    
                                    for (k=0;k<kernel_size;k++){
                                        
                                        beamsq_real[k] = 0.;
                                        beamsq_imag[k] = 0.;
                                        beam_real[k] = 0.;
                                        beamsq_real2[k] = 0.;
                                        beamsq_imag2[k] = 0.;
                                        beam_real2[k] = 0.;
                                        
                                        
                                        if (sqrt((u_lex_small[k]-uu)*(u_lex_small[k]-uu)+(v_lex_small[k]-vv)*(v_lex_small[k]-vv)) <= delt){
                                            
                                            for (jjj=0;jjj<4;jjj++){
                                                
                                                beamsq_real[k] += a_nut[jjj]*cos(2.*M_PI*((float) jjj)*(u_lex_small2[k]-uu2)/2./delt);
                                                beamsq_imag[k] += 0.;
                                                beam_real[k] += a_nut[jjj]*cos(2.*M_PI*((float) jjj)*(v_lex_small2[k]-vv2)/2./delt);
                                                beamsq_real2[k] += a_nut[jjj]*cos(2.*M_PI*((float) jjj)*(u_lex_small[k]-uu)/2./delt);
                                                beamsq_imag2[k] += 0.;
                                                beam_real2[k] += a_nut[jjj]*cos(2.*M_PI*((float) jjj)*(v_lex_small[k]-vv)/2./delt);
                                                
                                                
                                            }
                                            
                                            beamsq_real[k] *= beam_real[k];
                                            beamsq_real2[k] *= beam_real2[k];
                                            
                                            
                                                // new stuff for testing
                                            
                                            beamsq_real[k] = cabs(beamsq_real[k]);
                                            beamsq_real2[k] = cabs(beamsq_real2[k]);
                                            
                                                //
                                            
                                            
                                            interorig += beamsq_real[k];
                                                //        printf("Difference of beams %f: %f %f %f %f\n",beamsq_real[k]-beamsq_real2[k],uu,uu2,vv,vv2);
                                            
                                        }
                                        
                                        
                                    }
                                    
                                    
                                    
                                    total_norm_sq = creal(interorig);
                                    
                                        //total_norm_sq = 1.;
                                    
                                    normal = total_norm_sq;
                                    
                                        //printf("total_norm_sq: %f %f %f %f\n",creal(interorig),cimag(interorig),total_norm_sq,normal*DELTA_U*DELTA_U);
                                    
                                    
                                    
                                    /* Flag baselines with large amplitudes */
                                    
                                    
                                    
                                    vis_r1 = data->visdata[0][2*(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)];
                                    vis_r2 = data2->visdata[0][2*(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)];
                                    vis_i1 = data->visdata[0][2*(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1];
                                    vis_i2 = data2->visdata[0][2*(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1];
                                    
                                    
                                    flagbase = 0;
                                    
                                    if (sqrt((data->visdata[0][2*(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)])*(data->visdata[0][2*(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)]) + (data->visdata[0][2*(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1])*(data->visdata[0][2*(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1])) > 20000.){
                                        
                                        flagbase = 1;
                                        printf("Flagged baseline: %d\n",iii);
                                        
                                    }
                                    
                                    /* ONLY GRID THESE POINTS IF THE WEIGHTS MATCH */
                                    
                                    if (((data->weightdata[0][(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] == data2->weightdata[0][(iii*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)])&&(data->weightdata[0][(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] > 0.))&&(flagbase == 0)){
                                        
                                        
                                            //    printf("%f %f\n",vis_r1,vis_i1);
                                        
                                        /********************************************************/
                                            // Only process if finite
                                        if ((isnan(vis_r1+vis_r2+vis_i1+vis_i2+data->weightdata[0][(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)]) == 0)){
                                            /********************************************************/
                                            
                                                //        printf("%f %f %f %f %f %f %f\n",uu,vv,vis_r1,vis_i1,vis_r2,vis_i2,data->weightdata[0][(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)]);
                                            
                                                //if (sqrt(vis_r1*vis_r1+vis_i1*vis_i1) > 20000.) printf("i %duvw %f %f %f reals %lg %lg imag %lg %lg pol %d\n",i,uu,vv,ww,vis_r1,vis_r2,vis_i1,vis_i2,pol);
                                            /************************* TESTING ***********************/
                                            
                                            
                                            
                                            v1 = gauss_noise(1.)/sqrt(2.);
                                            v2 = gauss_noise(1.)/sqrt(2.);
                                            v3 = gauss_noise(1.)/sqrt(2.);
                                            v4 = gauss_noise(1.)/sqrt(2.);
                                                //           printf("v1 %f v2 %f\n",v1,v2);
                                            
                                            /* weights for each visibility */
                                            weight = (0.5*(data->weightdata[0][(iii*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] + data2->weightdata[0][(iii*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)]));
                                            
                                            if ((debug)&&(weight != 0.0)) fprintf(flog,"chan_width %f period %f weight %f reduced weight %f\n",global_chanwidth,global_period,weight,sqrt(weight/(global_chanwidth/1.e4*global_period/8.)));
                                            
                           /* Here we are trying to get the weights right */
                                            
                                            weight = sqrt(weight/64.);
                                            weight_norm = weight;
//                                            weight_norm=1;
                                    //        weight = 1.;
                                            
                                            
                                            
                                            /*********************************************/
                                            
                                                //    printf("weight %f\n",weight);
                                            
                                            
                                                //             printf("weight %f norms %f %f\n",weight,norm,normal);
                                            
                                                //            printf("pol %d abs vis %f\n",pol,vis_r*vis_r+vis_i*vis_i);
                                            
                                            if ((debug)&&(weight != 0.0)) fprintf(flog,"chan_width %f period %f ch %d pol %d weight %f vis_r %f beam %f loc %d\n",global_chanwidth,global_period,ch,pol,weight,creal(vis_itot),(beamsq_real[10]),xx);
                                            
                                            if (debug) fflush(flog);
                                            
                                            /* update bdag_b matrix and bdag_v vector according to the footprint of this baseline */
                                            
                                            
                                            
                                            
                                            /* START FIRST DATASET */
                                            
                                            
                                            for (jjj=0;jjj<kernel_size;jjj++){
                                                
                                                flag2 = 0;
                                                
                                                /* compute u and v locations for each cell */
                                                
                                                double u_loc1 = u_lex_small[jjj] -round(uu/DELTA_U)*DELTA_U + (uu);
                                                double v_loc1 = v_lex_small[jjj] -round(vv/DELTA_U)*DELTA_U + (vv);
                                                
                                                    //                  exit(1);
                                                
                                                /* determine index values for u location u_loc1 etc on THE FULL/HALF UV-PLANE */
                                                
                                                if (HALF_PLANE_FLAG) { ii = round(u_loc1/DELTA_U);} else { ii = round(u_loc1/DELTA_U) + u_size-1; }
                                                jj = round(v_loc1/DELTA_U) + u_size-1;
                                                
                                                /* If location goes into -ve u value, flag this as a location that needs to be re-gridded onto half-plane */
                                                if ((HALF_PLANE_FLAG)&&(ii < 0)){
                                                    ii = -ii;
                                                    jj = 2*u_size - jj;
                                                    flag2 = 1;
                                                }
                                                
                                                /* If location is within FULL uv-plane */
                                                if (HALF_PLANE_FLAG){ top = u_size; } else {top = 2.*u_size;}
                                                if ((ii >= 0)&&(ii <= top-1)&&(jj >= 0)&&(jj <= 2*u_size-1)){
                                                    
                                                    xx = round((ii)*2.*u_size + jj); /* lexicographic location */
                                                        //     if (xx<0) printf("xx ii jj %d %d %d\n",xx,ii,jj);
                                                    
                                                    
                                                    
                                                    if ((flag == 0)&&(flag2 == 0)){
                                                        
                                                        
                                                        bdag_v_real[xx][ch] += (vis_r1*beamsq_real2[jjj] )/normal/norm*weight*weight +I*( vis_i1*beamsq_real2[jjj] )/normal/norm*weight*weight;
                                                        
                                                        bdag_v_realdiff[xx][ch] += ( vis_r1*beamsq_real2[jjj] )/normal/norm*weight*weight +I*(vis_i1*beamsq_real2[jjj]  )/normal/norm*weight*weight;
                                                        
                                                        
                                                        
                                                    } else {
                                                        
                                                        
                                                        bdag_v_real[xx][ch] += (vis_r1*beamsq_real2[jjj])/normal/norm*weight*weight +I*( - vis_i1*beamsq_real2[jjj] )/normal/norm*weight*weight;
                                                        
                                                        bdag_v_realdiff[xx][ch] += ( vis_r1*beamsq_real2[jjj] )/normal/norm*weight*weight +I*(- vis_i1*beamsq_real2[jjj]  )/normal/norm*weight*weight;
                                                        
                                                        
                                                    }
                                                    if ((debug)&&(ch == 9)) fprintf(flog,"bdag %lf %lf weight %f\n",creal(bdag_v_real[xx][ch]),cimag(bdag_v_real[xx][ch]),weight);
                                                    
                                                    
                                                    /* EXPECTED NOISE CALC: Compute the expected contribution from a thermal-noise only visibility, located at the same point in the uvw-plane */
                                                    
                                                    
                                                    
                                                    if (flag2 == 0){
                                                        noise_array[xx][ch] += (beamsq_real2[jjj]*(v1+v3) + beamsq_imag2[jjj]*(v2+v4))/normal/norm*weight*weight + I*(-beamsq_imag2[jjj]*(v1+v3) + beamsq_real2[jjj]*(v2+v4))/normal/norm*weight*weight;
                                                            // printf("noise xx %f %d\n",noise_array[xx][ch],xx);
                                                    } else {
                                                        noise_array[xx][ch] += (beamsq_real2[jjj]*(v1+v3) + beamsq_imag2[jjj]*(v2+v4))/normal/norm*weight*weight - I*(-beamsq_imag2[jjj]*(v1+v3) + beamsq_real2[jjj]*(v2+v4))/normal/norm*weight*weight;
                                                    }
                                                    
                                                    if (flag2 == 0){
                                                        noisediff_array[xx][ch] += (beamsq_real2[jjj]*(v1-v3) + beamsq_imag2[jjj]*(v2-v4))/normal/norm*weight*weight + I*(-beamsq_imag2[jjj]*(v1-v3) + beamsq_real2[jjj]*(v2-v4))/normal/norm*weight*weight;
                                                    } else {
                                                        noisediff_array[xx][ch] += (beamsq_real2[jjj]*(v1-v3) + beamsq_imag2[jjj]*(v2-v4))/normal/norm*weight*weight - I*(-beamsq_imag2[jjj]*(v1-v3) + beamsq_real2[jjj]*(v2-v4))/normal/norm*weight*weight;
                                                    }
                                                    
                                                    
                                                    
                                                    
                                                    
                                                    flag_array[xx][ch] += ( (beamsq_real2[jjj]) )/normal/norm*weight_norm*weight_norm/2.;    //test_4
                                                    
                                                    
                                                    
                                                }
                                            }   /* end loop over beam */
                                            
                                            
                                                //        printf("End first dataset\n");
                                            
                                            /* END FIRST DATASET */
                                            
                                            
                                            /* SECOND DATASET */
                                            
                                            for (jjj=0;jjj<kernel_size;jjj++){
                                                
                                                flag2 = 0;
                                                
                                                /* compute u and v locations for each cell */
                                                
                                                double u_loc1 = u_lex_small2[jjj] -round(uu2/DELTA_U)*DELTA_U + (uu2);
                                                double v_loc1 = v_lex_small2[jjj] -round(vv2/DELTA_U)*DELTA_U + (vv2);
                                                
                                                    //                  exit(1);
                                                
                                                /* determine index values for u location u_loc1 etc on THE FULL/HALF UV-PLANE */
                                                
                                                if (HALF_PLANE_FLAG) { ii = round(u_loc1/DELTA_U);} else { ii = round(u_loc1/DELTA_U) + u_size-1; }
                                                jj = round(v_loc1/DELTA_U) + u_size-1;
                                                
                                                /* If location goes into -ve u value, flag this as a location that needs to be re-gridded onto half-plane */
                                                if ((HALF_PLANE_FLAG)&&(ii < 0)){
                                                    ii = -ii;
                                                    jj = 2*u_size - jj;
                                                    flag2 = 1;
                                                }
                                                
                                                /* If location is within FULL uv-plane */
                                                if (HALF_PLANE_FLAG){ top = u_size; } else {top = 2.*u_size;}
                                                if ((ii >= 0)&&(ii <= top-1)&&(jj >= 0)&&(jj <= 2*u_size-1)){
                                                    
                                                    xx = round((ii)*2.*u_size + jj); /* lexicographic location */
                                                    
                                                    if ((flag == 0)&&(flag2 == 0)){
                                                        
                                                        
                                                        bdag_v_real[xx][ch] += (vis_r2*beamsq_real[jjj] )/normal/norm*weight*weight +I*(+ vis_i2*beamsq_real[jjj] )/normal/norm*weight*weight;
                                                        
                                                        bdag_v_realdiff[xx][ch] += ( - vis_r2*beamsq_real[jjj] )/normal/norm*weight*weight +I*( - vis_i2*beamsq_real[jjj] )/normal/norm*weight*weight;
                                                        
                                                        
                                                    } else {
                                                        
                                                        bdag_v_real[xx][ch] += (vis_r2*beamsq_real[jjj] )/normal/norm*weight*weight +I*(  - vis_i2*beamsq_real[jjj] )/normal/norm*weight*weight;
                                                        
                                                        bdag_v_realdiff[xx][ch] += ( - vis_r2*beamsq_real[jjj] )/normal/norm*weight*weight +I*( + vis_i2*beamsq_real[jjj] )/normal/norm*weight*weight;
                                                        
                                                        
                                                    }
                                                    if ((debug)&&(ch == 9)) fprintf(flog,"bdag %lf %lf weight %f\n",creal(bdag_v_real[xx][ch]),cimag(bdag_v_real[xx][ch]),weight);
                                                    
                                                    
                                                    flag_array[xx][ch] += ( (beamsq_real[jjj]) )/normal/norm*weight_norm*weight_norm/2.;
                                                    
                                                    
                                                }
                                            }   /* end loop over beam */
                                            
                                            
                                                //           printf("End second dataset\n");
                                            
                                            /* END SECOND DATASET */
                                            
                                            
                                            
                                        } //end if statement for vis being finite
                                        
                                    }    /* end loop over weights being equal */
                                    
                                    
                                    
                                }   /* end IF/ELSE loop for baseline inclusion */
                                
                                
                                
                            } /* end loop over baselines for a given w stack and pol */
                            
                                //   if (addinternal > 0) printf("addinternal %lg %lg\n",additup,addinternal);
                            
                        } /* end if statement - don't process anything if no vis contribute!! */
                        
                        
                        
                        
                        /* write completion to log */
                        fprintf(flog,"Completed frequency index %d, channel %d, frequency %g, pol %d, at timestep %d and time %f with systemp %f\n",freq_index,ch,frequency,pol,chunk,data->date[0],systemp);
                            //      printf("Completed frequency index %d, channel %d, frequency %g, pol %d, at timestep %d and time %f with systemp %f\n",freq_index,ch,frequency,pol,chunk,data->date[0],systemp);
                            //fflush(flog);
                        
                        
                        /* clean-up */
                        free(beam_real);
                        free(beam_imag);
                        free(beam_real2);
                        free(beam_imag2);
                        free(beamsq_real);
                        free(beamsq_imag);
                        free(beamsq_real2);
                        free(beamsq_imag2);
                        free(u_lex_small);
                        free(v_lex_small);
                        free(u_lex_small2);
                        free(v_lex_small2);
                        free(include_vis);
                        
                        fflush(flog);
                        
                    }  //  ******** end loop over channels ********
                    
                    
                    
                    
                        //             printf("Freed \n");
                    
                    
                        //       printf("Completed timestep %d and time %f with systemp %f and channels %d\n",chunk,data->date[0],systemp,data->n_freq);
                    
                }  /*************** END OPENMP PARALLEL LOOP ********************/
                
                
                
                
                
                if ((res=readUVFITSnextIter(data,iter)) !=0){
                    printf("First set iteration not found at the end\n");
                    res=0;
                    ffcmrk;
                    break;
                }
                
                if ((res=readUVFITSnextIter(data2,iter2)) !=0){
                    printf("Second set iteration not found at the end\n");
                    res=0;
                    ffcmrk;
                    break;
                }
                
                
            } else {
                
                printf("Second set does not exist, therefore not gridding final odd timestep\n");
                
            }// end if statement to only run if the second set is present
            
        }  /*************** END WHILE LOOP OVER UVFITS TIME STEPS *******/
        /*       if(res != 1) {
         fprintf(stderr,"readUVFISTnextIter returned %d\n",res);
         }
         */
        
        printf("Do I make it here finally?\n");
        /* write-out output files */
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        /************** OUTPUT FILES ********************/
        /* define filenames */
        
        JD_to_Cal(data->date[0],&obs_year, &obs_month, &obs_day);
            //    printf("Year %d Month %d Day %d\n",obs_year,obs_month,obs_day);
        
        if (pol == 0){
            sprintf(filename_real,"%snoisec_xx.%s.dat",getenv("OUTPUTDIR"),ext);
            sprintf(filename_real2,"%sbv_xx.%s.dat",getenv("OUTPUTDIR"),ext);
            sprintf(filename_reald,"%snoisecdiff_xx.%s.dat",getenv("OUTPUTDIR"),ext);
            sprintf(filename_real2d,"%sbvdiff_xx.%s.dat",getenv("OUTPUTDIR"),ext);
            sprintf(filename_flags,"%sweightsc_xx.%s.dat",getenv("OUTPUTDIR"),ext);
        } else {
            sprintf(filename_real,"%snoisec_yy.%s.dat",getenv("OUTPUTDIR"),ext);
            sprintf(filename_real2,"%sbv_yy.%s.dat",getenv("OUTPUTDIR"),ext);
            sprintf(filename_reald,"%snoisecdiff_yy.%s.dat",getenv("OUTPUTDIR"),ext);
            sprintf(filename_real2d,"%sbvdiff_yy.%s.dat",getenv("OUTPUTDIR"),ext);
            sprintf(filename_flags,"%sweightsc_yy.%s.dat",getenv("OUTPUTDIR"),ext);
        }
        
        
        printf("Writing output files...\n");
        
        
        
            // Check to see if the file already exists
        
        if ((fptr = fopen(filename_real,"r")) == NULL){
            
            if ((fptr = fopen(filename_real,"w")) != NULL){
                for (i=0;i<size_fin;i++){
                    for (j=0;j<data->n_freq;j++){
                        fwrite(&noise_array[i][j],sizeof(double complex),1,fptr);
                    }
                }
                fclose(fptr);
            } else {
                fprintf(flog,"Failed to write observation noise file: %s\n",filename_real);
            }
            
        } else {
            
            printf("File exists... opening...\n");
            
                // temp_noise = malloc(total_size*sizeof(double complex));
            temp_noise = create2DMatrixComplex(size_fin,data->n_freq);
            assert(temp_noise != NULL);
            for (i=0;i<size_fin;i++){
                for (j=0;j<data->n_freq;j++){
                    fread(&temp_noise[i][j],sizeof(double complex),1,fptr);
                }
            }
            
            fclose(fptr);
            
            for (i=0;i<size_fin;i++){
                for (j=0;j<data->n_freq;j++){
                    noise_array[i][j] += temp_noise[i][j];
                }
            }
            
            free_memComplex(temp_noise,size_fin,data->n_freq);
            
                // write_out data
            
            if ((fptr = fopen(filename_real,"w")) != NULL){
                for (i=0;i<size_fin;i++){
                    for (j=0;j<data->n_freq;j++){
                        fwrite(&noise_array[i][j],sizeof(double complex),1,fptr);
                    }
                }
                
                
            } else {
                fprintf(flog,"Failed to increment observation noise file: %s\n",filename_real);
            }
            
        }
        
        
        
        free_memComplex(noise_array,size_fin,data->n_freq);
        
        
            // Check to see if the file already exists
        
        if ((fptr = fopen(filename_real2,"r")) == NULL){
            
            if ((fptr = fopen(filename_real2,"w")) != NULL){
                for (i=0;i<size_fin;i++){
                    for (j=0;j<data->n_freq;j++){
                        fwrite(&bdag_v_real[i][j],sizeof(double complex),1,fptr);
                    }
                }
                fclose(fptr);
            } else {
                fprintf(flog,"Failed to write observation vis tot file: %s\n",filename_real2);
            }
            
        } else {
            
            printf("File exists... opening...\n");
            
            temp_v_real = create2DMatrixComplex(size_fin,data->n_freq);
            assert(temp_v_real != NULL);
            for (i=0;i<size_fin;i++){
                for (j=0;j<data->n_freq;j++){
                    fread(&temp_v_real[i][j],sizeof(double complex),1,fptr);
                }
            }
            
            
            
            fclose(fptr);
            
            for (i=0;i<size_fin;i++){
                for (j=0;j<data->n_freq;j++){
                    bdag_v_real[i][j] += temp_v_real[i][j];
                }
            }
            
            free_memComplex(temp_v_real,size_fin,data->n_freq);
            
                //     free(temp_noise);
            
                // write_out data
            
            if ((fptr = fopen(filename_real2,"w")) != NULL){
                for (i=0;i<size_fin;i++){
                    for (j=0;j<data->n_freq;j++){
                        fwrite(&bdag_v_real[i][j],sizeof(double complex),1,fptr);
                    }
                }
                
            } else {
                fprintf(flog,"Failed to increment observation noise file: %s\n",filename_real);
            }
            
            
        }
        
        
        free_memComplex(bdag_v_real,size_fin,data->n_freq);
        
        /* DIFF FILES */
        
        
            // Check to see if the file already exists
        
        if ((fptr = fopen(filename_reald,"r")) == NULL){
            
            if ((fptr = fopen(filename_reald,"w")) != NULL){
                for (i=0;i<size_fin;i++){
                    for (j=0;j<data->n_freq;j++){
                        fwrite(&noisediff_array[i][j],sizeof(double complex),1,fptr);
                    }
                }
                fclose(fptr);
                
            } else {
                fprintf(flog,"Failed to write observation noise file: %s\n",filename_reald);
            }
            
        } else {
            
            printf("File exists... opening...\n");
            
            temp_noisediff = create2DMatrixComplex(size_fin,data->n_freq);
            assert(temp_noisediff != NULL);
            for (i=0;i<size_fin;i++){
                for (j=0;j<data->n_freq;j++){
                    fread(&temp_noisediff[i][j],sizeof(double complex),1,fptr);
                    
                }
            }
            
            
            fclose(fptr);
            
            
            
            for (i=0;i<size_fin;i++){
                for (j=0;j<data->n_freq;j++){
                    noisediff_array[i][j] += temp_noisediff[i][j];
                }
            }
            
            free_memComplex(temp_noisediff,size_fin,data->n_freq);
            
                //     free(temp_noise);
            
                // write_out data
            
            if ((fptr = fopen(filename_reald,"w")) != NULL){
                for (i=0;i<size_fin;i++){
                    for (j=0;j<data->n_freq;j++){
                        fwrite(&noisediff_array[i][j],sizeof(double complex),1,fptr);
                    }
                }
                
            } else {
                fprintf(flog,"Failed to increment observation noise file: %s\n",filename_reald);
            }
            
            
        }
        
        
        
        free_memComplex(noisediff_array,size_fin,data->n_freq);
        
        
        
            // Check to see if the file already exists
        
        if ((fptr = fopen(filename_real2d,"r")) == NULL){
            
            if ((fptr = fopen(filename_real2d,"w")) != NULL){
                for (i=0;i<size_fin;i++){
                    for (j=0;j<data->n_freq;j++){
                        fwrite(&bdag_v_realdiff[i][j],sizeof(double complex),1,fptr);
                    }
                }
                fclose(fptr);
            } else {
                fprintf(flog,"Failed to write observation vis tot file: %s\n",filename_real2d);
            }
            
        } else {
            
                //        printf("File exists... opening...\n");
            
            temp_v_realdiff = create2DMatrixComplex(size_fin,data->n_freq);
            assert(temp_v_realdiff != NULL);
            for (i=0;i<size_fin;i++){
                for (j=0;j<data->n_freq;j++){
                    fread(&temp_v_realdiff[i][j],sizeof(double complex),1,fptr);
                }
            }
            
            
            
            fclose(fptr);
            
            
            
            for (i=0;i<size_fin;i++){
                for (j=0;j<data->n_freq;j++){
                    bdag_v_realdiff[i][j] += temp_v_realdiff[i][j];
                }
            }
            
            free_memComplex(temp_v_realdiff,size_fin,data->n_freq);
            
                //     free(temp_noise);
            
                // write_out data
            
            if ((fptr = fopen(filename_real2d,"w")) != NULL){
                for (i=0;i<size_fin;i++){
                    for (j=0;j<data->n_freq;j++){
                        fwrite(&bdag_v_realdiff[i][j],sizeof(double complex),1,fptr);
                    }
                }
                
                
            } else {
                fprintf(flog,"Failed to increment observation noise file: %s\n",filename_real2d);
            }
            
            
        }
        
        
        
        
        free_memComplex(bdag_v_realdiff,size_fin,data->n_freq);
        
        
        /* WEIGHTS FILE */
        
        /* write-out flags file of sampled uv locations */
        
        if ((fptrflags = fopen(filename_flags,"r")) == NULL){
            
            if ((fptrflags = fopen(filename_flags,"w")) != NULL){
                for (i=0;i<size_fin;i++){
                    for (j=0;j<data->n_freq;j++){
                        fwrite(&flag_array[i][j],sizeof(double complex),1,fptrflags);
                    }
                }
                
                fclose(fptrflags);
            } else {
                fprintf(flog,"Failed to write observation flags file: %s\n",filename_flags);
            }
            
        } else {
            
            temp_flags = create2DMatrixComplex(size_fin,data->n_freq);
            assert(temp_flags != NULL);
            for (i=0;i<size_fin;i++){
                for (j=0;j<data->n_freq;j++){
                    fread(&temp_flags[i][j],sizeof(double complex),1,fptrflags);
                    
                }
            }
            
            
            fclose(fptrflags);
            
            
            
            for (i=0;i<size_fin;i++){
                for (j=0;j<data->n_freq;j++){
                    flag_array[i][j] += temp_flags[i][j];
                }
            }
            
            
                //  free(temp_flags);
            
                // write_out data
            
            if ((fptrflags = fopen(filename_flags,"w")) != NULL){
                for (i=0;i<size_fin;i++){
                    for (j=0;j<data->n_freq;j++){
                        fwrite(&flag_array[i][j],sizeof(double complex),1,fptrflags);
                    }
                }
            } else {
                fprintf(flog,"Failed to increment observation flags file: %s\n",filename_flags);
            }
            
            
            
        }
        
        free_memComplex(flag_array,size_fin,data->n_freq);
        
        readUVFITSCloseIter(iter);
        readUVFITSCloseIter(iter2);
        
        
    }  // end loop over polarisations
        
    fclose(flog);

    
 
    
    // end if statement for band selection
    

    return 0;
}





// Functions


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


    double complex **create2DMatrixComplex(int sizex, int sizey){
    double complex **output;
    int k;

        output = calloc(sizex,sizeof(double complex*));
        for (k=0;k<sizex;k++){
            output[k] = calloc(sizey,sizeof(double complex));
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


void memzero3Dc(double complex ***matrix,int xsize, int ysize, int zsize){
  int i,j;

  for (i=0;i<xsize;i++){
    for (j=0;j<ysize;j++) memset(matrix[i][j],0.,zsize*sizeof(double complex));
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

void free_memComplex(double complex **matrix, int xsize, int ysize){
    int i;
    
    for (i=0;i<xsize;i++) free(matrix[i]);
    
}


void free_mem3D(double ***matrix, int xsize, int ysize, int zsize){
    int i,j;
    
    for (i=0;i<xsize;i++){
        for (j=0;j<ysize;j++) free(matrix[i][j]);
            }

}




