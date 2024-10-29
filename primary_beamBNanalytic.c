/************************************
 ! MODULE: primary_beam
 !         utilities for reading MWA primary beam data
 !         Author: Cathryn Trott. August 2014
 !		Version: 1.3
 ************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <float.h>
#include <assert.h>
#include "uvfits.h"
#include "assert.h"
#include "primary_beam.h"
#include "fitsio.h"
#include "complex.h"

/*
extern float DELTA_WSTACK;
extern float DELTA_W;
extern float DELTA_U;
extern float UMAX;
extern float PERIOD;
*/

/* private function prototypes */
int beam_diff_mwa_vector_2D(float **beam_r, float **beam_i, float *u_tab, double *u,double uprime,double *v,double vprime,int size,double *out_real,double *out_imag, float factor_scale);
int interp(float **uvgrid, double u, double v, double *res_re);
int interp_bilin(float **uvgrid, double u, double v, double *res_re);
double gauss_noise (double value1);
double modd(double val);

/* private global vars */
static int debug=0;

/********************************************************************************/

/******************************
 ! NAME:		calc_bv_flags
 ! PURPOSE:		compute the bdag_b and bdag_v output files for each freq, pol and w-stack. centres beam on uv-cell. Small w-stacks
 ! ARGUMENTS:	
 ! RETURNS:		integer (0=success)
******************************/

int calc_bv_bn_analytic(uvdata *data, uvdata *data2, int size_fin, int ch, double frequency, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext, int *flag14m, int band, float umx){

FILE *fptr=NULL,*fptrv=NULL,*fptrflags=NULL,*flog=NULL;
int flag=0,flag2=0,w_index,obs_year,obs_month,obs_day,i,j,k,ii,jj,xx,beam_size,kernel_size,top,*include_vis=NULL,loc=0,flagbase=0;
double w_centre,checksum_real=0.0,checksum_imag=0.0,checksum_abs=0.0,checksum_cent=0.0,v1,v2,v3,v4,beam_sigma2=0.05;
long contrib_vis=0,total_vis=0;
float vis_r1,vis_r2,vis_i1,vis_i2,rotat1,rotat2,rotat3,delt=DELTAUU;
char filename_real[1024],filename_real2[1024],filename_reald[1024],filename_real2d[1024],filename_flags[1024],debugfile[1024];
double complex  *bdag_v_real=NULL,*bdag_v_realdiff=NULL,*temp_v_real=NULL,*temp_noise=NULL,*temp_v_realdiff=NULL,*temp_noisediff=NULL,*noise_array=NULL,*noisediff_array=NULL,*flag_array=NULL,*temp_flags=NULL, inter=0.+I*0.,interorig=0.+I*0.;
    float vis_rtot,vis_rdiff,vis_itot,vis_idiff,w_temp=0.,kernel_size_float,visnoise,weight,normal=0.,total_norm_sq=0.,frequency_lower=0.,factor_scale=1.,l8,wetot=0.,weall=0.,additup=0.;
double distance=0.,*beam_real=NULL,*beam_imag=NULL,*u_lex_small=NULL,*v_lex_small=NULL,uu,vv,ww,norm,*beamsq_real=NULL,*beamsq_imag=NULL,*beamsq_real2=NULL,*beamsq_imag2=NULL,addinternal=0.,uu2,vv2,ww2,*u_lex_small2=NULL,*v_lex_small2=NULL,*beam_real2=NULL,*beam_imag2=NULL;
    double a_nut[4]={0.3635819,0.4891775,0.1365995,0.0106411};
    float nut_norm = 0.;
    
    
    
    a_nut[0]=0.3635819;
    a_nut[1]=0.4891775;
    a_nut[2]=0.1365995;
    a_nut[3]=0.0106411;   // B-N
    
       a_nut[0]=0.36358192677076108;
	a_nut[1]=0.489177437145017;
	a_nut[2]=0.136599513978692;
	a_nut[3]=0.01064112210553003;   // B-N
    


    srand48((long)time(NULL));
    
if (debug){
/* Open log file */
sprintf(debugfile,"%ssysdebug.txt",getenv("OUTPUTDIR"));
if ((flog=fopen(debugfile,"a")) == NULL){
	fprintf(stderr,"Cannot open output debug log file\n");
	return 1;
}
}

float intrinsic_res = 1./6.;

visnoise = 1.;
norm=1.;
    
/* define frequency at low end of coarse channel */
    
    frequency_lower = freq_index*COARSE_CHAN_WIDTH + LOWER_FREQ;
    factor_scale = frequency/frequency_lower;
//    factor_scale = 1.;
  //  printf("factor_scale %f\n",factor_scale);

/* Define dense vector for bdag_v accumulation and flags for uv-sampling */

	bdag_v_real = calloc(size_fin,sizeof(double complex));
	noise_array = calloc(size_fin,sizeof(double complex));
    bdag_v_realdiff = calloc(size_fin,sizeof(double complex));
    noisediff_array = calloc(size_fin,sizeof(double complex));
	flag_array = calloc(size_fin,sizeof(double complex));

for (wstack=0;wstack<NUM_W_STACK;wstack++){
    
	w_centre = (wstack*DELTA_WSTACK);
    
//	printf("wstack %d w_centre %lg\n",wstack,w_centre);

   /*****************************************************************/
   /* Loop over baselines - first loop to compute memory allocation */

	include_vis = calloc(data->n_baselines[0],sizeof(int));

	for (i=0;i<data->n_baselines[0];i++){
       //  printf("baseline: %d\n",i);
				
     /* Compute whether uv point lies within range of uv-grid and wstack */

		distance = sqrt(data->u[0][i]*data->u[0][i]+data->v[0][i]*data->v[0][i])*(frequency);
	//	if (debug) fprintf(flog,"i %d distance: %lg w: %g, weight %g, WMAX %f, wcentre %f\n",i,distance,data->w[0][i]*frequency,data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)],WMAX,w_centre);
    //    printf("Distance: %f\n",distance);

        
 		if ((distance >= 0.95*(umx-SIZE_BEAM/2*INTRINSIC_DELTA_U)) || (distance < 0.5) || (sqrt(data->w[0][i]*frequency*data->w[0][i]*frequency) >= WMAX) || ((sqrt((sqrt(data->w[0][i]*data->w[0][i])*frequency-w_centre)*(sqrt(data->w[0][i]*data->w[0][i])*frequency-w_centre))) > WSTACK_DELT_TOT/2.) || (data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] <= 0.) || (flag14m[i] == 0)){


		} else {

			contrib_vis++;
			include_vis[i] = 1;

		}




	}



	if (debug) fprintf(flog,"contrib vis: %ld\n",contrib_vis);

	/* set-up beam to contain enough entries to comfortably contain beam
			  region */
	
    beam_size = (int) BEAM_SIZE_FLOAT;
  //  printf("beam size %d\n",beam_size);
		
	if (((float)beam_size)/2. == (int)((float)(beam_size)/2.)) beam_size = beam_size+1;  /* make odd */

	kernel_size_float = beam_size*beam_size;
	kernel_size = (int)kernel_size_float;


   /*****************************************************************/
   /* Loop over baselines - second loop to actually process data */

	if (contrib_vis == 0){
	} else {

  addinternal=0.;
additup=0.;

	for (i=0;i<data->n_baselines[0];i++){
     //    printf("baseline: %d\n",i);
				

		if (include_vis[i] != 1){
		} else {

			//contrib_vis++;

			uu = data->u[0][i]*(frequency);
			vv = data->v[0][i]*(frequency);
			ww = data->w[0][i]*(frequency);
            uu2 = data2->u[0][i]*(frequency);
            vv2 = data2->v[0][i]*(frequency);
            ww2 = data2->w[0][i]*(frequency);
			flag = 0;
			

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
            
            u_lex_small = calloc(kernel_size,sizeof(double));
            v_lex_small = calloc(kernel_size,sizeof(double));
            u_lex_small2 = calloc(kernel_size,sizeof(double));
            v_lex_small2 = calloc(kernel_size,sizeof(double));

            loc = 0;
            
    			   /* lexicographic u and v co-ordinates of small *sky* Fourier plane around the Fourier location */
            for (k=0;k<beam_size;k++){
                for (j=0;j<beam_size;j++){
                    
                    u_lex_small[loc] = (double) ((k)-(int) beam_size/2)*DELTA_U + round(uu/DELTA_U)*DELTA_U;
                    v_lex_small[loc] = (double) ((j)-(int) beam_size/2)*DELTA_U + round(vv/DELTA_U)*DELTA_U;
                    u_lex_small2[loc] = (double) ((k)-(int) beam_size/2)*DELTA_U + round(uu2/DELTA_U)*DELTA_U;
                    v_lex_small2[loc] = (double) ((j)-(int) beam_size/2)*DELTA_U + round(vv2/DELTA_U)*DELTA_U;
                    if (debug) printf("loc %ld u_lex %g v_lex %g\n",loc,u_lex_small[loc],v_lex_small[loc]);
                    loc++;
                }
            }
            
			beam_real = calloc(kernel_size,sizeof(double));
			beam_imag = calloc(kernel_size,sizeof(double));

            beam_real2 = calloc(kernel_size,sizeof(double));
            beam_imag2 = calloc(kernel_size,sizeof(double));
            
			beamsq_real = calloc(kernel_size,sizeof(double));
			beamsq_imag = calloc(kernel_size,sizeof(double));
            
            beamsq_real2 = calloc(kernel_size,sizeof(double));
            beamsq_imag2 = calloc(kernel_size,sizeof(double));

			/* GET PRIMARY BEAM VALUES HERE */

            // use the values when the beam has not been squared. We are going to grid the vis with phased beam, and also grid the weights with the same complex phased beam. The gridded vis (not the weights) are then normalised by the square of the total of this phased beam to restore the full flux density of the sky.
            
//            get_pbsq_values_fine(beam_data,wstack,point,kernel_size,u_lex_small,v_lex_small,uu,vv,beamsq_real,beamsq_imag,1.,ch);


/* Compute beam gridding by an analytic Blackman-Harris, independent of frequency, w and pointing */

//	DELTAUU = 5.;
            if (band == 1)  delt = DELTAUU;
            if (band == 0)  delt = DELTAUU;//*152./182.;
            if ((band != 0)&&(band != 1))  delt = DELTAUU*82./182.;
          interorig = 0. + I*0.;
	
	for (k=0;k<kernel_size;k++){
        
        if (sqrt((u_lex_small[k]-uu)*(u_lex_small[k]-uu)+(v_lex_small[k]-vv)*(v_lex_small[k]-vv)) <= delt){
        
        for (j=0;j<4;j++){

		beamsq_real[k] += a_nut[j]*cos(2.*M_PI*((float) j)*(u_lex_small2[k]-uu2)/2./delt);
		beamsq_imag[k] += 0.;
            beam_real[k] += a_nut[j]*cos(2.*M_PI*((float) j)*(v_lex_small2[k]-vv2)/2./delt);
		beamsq_real2[k] += a_nut[j]*cos(2.*M_PI*((float) j)*(u_lex_small[k]-uu)/2./delt);
		beamsq_imag2[k] += 0.;
            beam_real2[k] += a_nut[j]*cos(2.*M_PI*((float) j)*(v_lex_small[k]-vv)/2./delt);
        
            
        }
        
        beamsq_real[k] *= beam_real[k];
        beamsq_real2[k] *= beam_real2[k];
        
            
            // new stuff for testing
            
            beamsq_real[k] = cabs(beamsq_real[k]);
            beamsq_real2[k] = cabs(beamsq_real2[k]);
            
            //
  
            
            interorig += beamsq_real[k];
        //    printf("Difference of beams %f: %f %f %f %f\n",beamsq_real[k]-beamsq_real2[k],uu,uu2,vv,vv2);
        
        }
  

	}


            
            total_norm_sq = creal(interorig);

//total_norm_sq = 1.;
            
            normal = total_norm_sq;
            
  //          printf("total_norm_sq: %f %f %f %f\n",creal(interorig),cimag(interorig),total_norm_sq,normal*DELTA_U*DELTA_U);
            
                 
            
            /* Flag baselines with large amplitudes */
            
            flagbase = 0;
            
            if (sqrt((data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)])*(data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)]) + (data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1])*(data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1])) > 20000.){
                
            flagbase = 1;
                printf("Flagged baseline: %d\n",i);
                
            }
            
            
        /* ONLY GRID THESE POINTS IF THE WEIGHTS MATCH */
            
           if (((data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] == data2->weightdata[0][(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)])&&(data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] > 0.))&&(flagbase == 0)){
            
     
            
            vis_r1 = data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)];
            vis_r2 = data2->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)];
            vis_i1 = data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1];
            vis_i2 = data2->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1];


	//   printf("%f %f\n",vis_r1,vis_i1);

/********************************************************/
// Only process if finite			   
if ((isnan(vis_r1+vis_r2+vis_i1+vis_i2+data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)]) == 0)){
/********************************************************/
			  
//		printf("%f %f %f %f %f %f %f\n",uu,vv,vis_r1,vis_i1,vis_r2,vis_i2,data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)]);

if (sqrt(vis_r1*vis_r1+vis_i1*vis_i1) > 20000.) printf("i %duvw %f %f %f reals %lg %lg imag %lg %lg pol %d\n",i,uu,vv,ww,vis_r1,vis_r2,vis_i1,vis_i2,pol);
/************************* TESTING ***********************/



            v1 = gauss_noise(1.)/sqrt(2.);
            v2 = gauss_noise(1.)/sqrt(2.);
            v3 = gauss_noise(1.)/sqrt(2.);
            v4 = gauss_noise(1.)/sqrt(2.);
 //           printf("v1 %f v2 %f\n",v1,v2);

			/* weights for each visibility */
			weight = 0.5*(data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] + data2->weightdata[0][(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)]);

		if ((debug)&&(weight != 0.0)) fprintf(flog,"chan_width %f period %f weight %f reduced weight %f\n",global_chanwidth,global_period,weight,sqrt(weight/(global_chanwidth/1.e4*global_period/8.)));


            weight = sqrt(weight/(64.));


//	printf("weight %f\n",weight);


               wetot += weight;
               weall += 1.;
               
            //  printf("weight %f norms %f %f\n",weight,norm,normal);
               
//            printf("pol %d abs vis %f\n",pol,vis_r*vis_r+vis_i*vis_i);

            if ((debug)&&(weight != 0.0)) fprintf(flog,"chan_width %f period %f ch %d pol %d weight %f vis_r %f beam %f loc %d\n",global_chanwidth,global_period,ch,pol,weight,creal(vis_itot),(beamsq_real[10]),xx);

			if (debug) fflush(flog);

		       /* update bdag_b matrix and bdag_v vector according to the footprint of this baseline */
               
               
               
               
               /* START FIRST DATASET */
             

               for (j=0;j<kernel_size;j++){
                   
                   flag2 = 0;
                   
                   /* compute u and v locations for each cell */
                   
                   double u_loc1 = u_lex_small[j] -round(uu/DELTA_U)*DELTA_U + (uu);
                   double v_loc1 = v_lex_small[j] -round(vv/DELTA_U)*DELTA_U + (vv);
                   
                   //      			exit(1);
                   
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
                           
                          
                           bdag_v_real[xx] += (vis_r1*beamsq_real2[j] )/normal/norm*weight*weight +I*( vis_i1*beamsq_real2[j] )/normal/norm*weight*weight;
                           
                           bdag_v_realdiff[xx] += ( vis_r1*beamsq_real2[j] )/normal/norm*weight*weight +I*(vis_i1*beamsq_real2[j]  )/normal/norm*weight*weight;
                           
                           
                           
                       } else {
                           
                          
                           bdag_v_real[xx] += (vis_r1*beamsq_real2[j])/normal/norm*weight*weight +I*( - vis_i1*beamsq_real2[j] )/normal/norm*weight*weight;
                           
                           bdag_v_realdiff[xx] += ( vis_r1*beamsq_real2[j] )/normal/norm*weight*weight +I*(- vis_i1*beamsq_real2[j]  )/normal/norm*weight*weight;
                           
                           
                       }
                       if ((debug)&&(ch == 9)) fprintf(flog,"bdag %lf %lf weight %f\n",creal(bdag_v_real[xx]),cimag(bdag_v_real[xx]),weight);
                       
                       
                       /* EXPECTED NOISE CALC: Compute the expected contribution from a thermal-noise only visibility, located at the same point in the uvw-plane */
                       
                       
                       
                       if (flag2 == 0){
                           noise_array[xx] += (beamsq_real2[j]*(v1+v3) + beamsq_imag2[j]*(v2+v4))/normal/norm*weight*weight + I*(-beamsq_imag2[j]*(v1+v3) + beamsq_real2[j]*(v2+v4))/normal/norm*weight*weight;
                          // printf("noise xx %f %d\n",noise_array[xx],xx);
                       } else {
                           noise_array[xx] += (beamsq_real2[j]*(v1+v3) + beamsq_imag2[j]*(v2+v4))/normal/norm*weight*weight - I*(-beamsq_imag2[j]*(v1+v3) + beamsq_real2[j]*(v2+v4))/normal/norm*weight*weight;
                       }
                       
                       if (flag2 == 0){
                           noisediff_array[xx] += (beamsq_real2[j]*(v1-v3) + beamsq_imag2[j]*(v2-v4))/normal/norm*weight*weight + I*(-beamsq_imag2[j]*(v1-v3) + beamsq_real2[j]*(v2-v4))/normal/norm*weight*weight;
                       } else {
                           noisediff_array[xx] += (beamsq_real2[j]*(v1-v3) + beamsq_imag2[j]*(v2-v4))/normal/norm*weight*weight - I*(-beamsq_imag2[j]*(v1-v3) + beamsq_real2[j]*(v2-v4))/normal/norm*weight*weight;
                       }
                       
                       
                       
 
                       
                       flag_array[xx] += ( (beamsq_real2[j]) )/normal/norm*weight*weight/2.;    //test_4
                       
            
                       
                   }
               }   /* end loop over beam */
               
               
       //        printf("End first dataset\n");
               
               /* END FIRST DATASET */
               

               /* SECOND DATASET */
               
               for (j=0;j<kernel_size;j++){
                   
                   flag2 = 0;
                   
                   /* compute u and v locations for each cell */
                   
                   double u_loc1 = u_lex_small2[j] -round(uu2/DELTA_U)*DELTA_U + (uu2);
                   double v_loc1 = v_lex_small2[j] -round(vv2/DELTA_U)*DELTA_U + (vv2);
                   
                   //      			exit(1);
                   
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
                           
                           
                           bdag_v_real[xx] += (vis_r2*beamsq_real[j] )/normal/norm*weight*weight +I*(+ vis_i2*beamsq_real[j] )/normal/norm*weight*weight;
                           
                           bdag_v_realdiff[xx] += ( - vis_r2*beamsq_real[j] )/normal/norm*weight*weight +I*( - vis_i2*beamsq_real[j] )/normal/norm*weight*weight;
                           
                           
                       } else {
                           
                           bdag_v_real[xx] += (vis_r2*beamsq_real[j] )/normal/norm*weight*weight +I*(  - vis_i2*beamsq_real[j] )/normal/norm*weight*weight;
                           
                           bdag_v_realdiff[xx] += ( - vis_r2*beamsq_real[j] )/normal/norm*weight*weight +I*( + vis_i2*beamsq_real[j] )/normal/norm*weight*weight;
                           
                           
                       }
                       if ((debug)&&(ch == 9)) fprintf(flog,"bdag %lf %lf weight %f\n",creal(bdag_v_real[xx]),cimag(bdag_v_real[xx]),weight);
                       
                       
                       flag_array[xx] += ( (beamsq_real[j]) )/normal/norm*weight*weight/2.;
                      
                       
                   }
               }   /* end loop over beam */
               
               
    //           printf("End second dataset\n");
               
               /* END SECOND DATASET */



            } //end if statement for vis being finite   
                
                
                
                
                
	}    /* end loop over weights being equal */
       

			total_vis++;
			if (debug) fprintf(flog,"total vis %ld\n",total_vis);

            
    //        printf("Freeing \n");
        
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

     //      printf("Freed \n");

		}   /* end IF/ELSE loop for baseline inclusion */



	} /* end loop over baselines for a given w stack and pol */
        
  //   if (addinternal > 0) printf("addinternal %lg %lg\n",additup,addinternal);

	} /* end if statement - don't process anything if no vis contribute!! */

	free(include_vis);

//	printf("Number of baselines contributing %ld of %d, at w-centre %g, polarization %d and frequency %g\n",contrib_vis,data->n_baselines[0],w_centre,pol,frequency);

    } /* End loop over wstacks before writing out files */
    
	if (debug) fclose(flog);
//	if (debug) exit(1);
    
 //   printf("Included weights %f of %f with %f\n",wetot,weall,additup);

	if (contrib_vis > 0){


/************** OUTPUT FILES ********************/
	/* define filenames */

	JD_to_Cal(data->date[0],&obs_year, &obs_month, &obs_day);
//	printf("Year %d Month %d Day %d\n",obs_year,obs_month,obs_day);

	if (pol == 0){
	sprintf(filename_real,"%snoise_freq%03.7g_xx.%s.dat",getenv("OUTPUTDIR"),frequency/1.e6,ext);
	sprintf(filename_real2,"%sbv_freq%03.7g_xx.%s.dat",getenv("OUTPUTDIR"),frequency/1.e6,ext);
    sprintf(filename_reald,"%snoisediff_freq%03.7g_xx.%s.dat",getenv("OUTPUTDIR"),frequency/1.e6,ext);
    sprintf(filename_real2d,"%sbvdiff_freq%03.7g_xx.%s.dat",getenv("OUTPUTDIR"),frequency/1.e6,ext);
	sprintf(filename_flags,"%sweights_freq%03.7g_xx.%s.dat",getenv("OUTPUTDIR"),frequency/1.e6,ext);
	} else {
	sprintf(filename_real,"%snoise_freq%03.7g_yy.%s.dat",getenv("OUTPUTDIR"),frequency/1.e6,ext);
	sprintf(filename_real2,"%sbv_freq%03.7g_yy.%s.dat",getenv("OUTPUTDIR"),frequency/1.e6,ext);
    sprintf(filename_reald,"%snoisediff_freq%03.7g_yy.%s.dat",getenv("OUTPUTDIR"),frequency/1.e6,ext);
    sprintf(filename_real2d,"%sbvdiff_freq%03.7g_yy.%s.dat",getenv("OUTPUTDIR"),frequency/1.e6,ext);
	sprintf(filename_flags,"%sweights_freq%03.7g_yy.%s.dat",getenv("OUTPUTDIR"),frequency/1.e6,ext);
	}

	// Check to see if the file already exists

	if ((fptr = fopen(filename_real,"r")) == NULL){

		if ((fptr = fopen(filename_real,"w")) != NULL){
			fwrite(noise_array,sizeof(double complex),size_fin,fptr);
			fclose(fptr);
		} else {
		fprintf(flog,"Failed to write observation noise file: %s\n",filename_real);
		}

	} else {

//		printf("File exists... opening...\n");

		temp_noise = malloc(size_fin*sizeof(double complex));
        assert(temp_noise != NULL);
        if (fread(temp_noise,sizeof(double complex),size_fin,fptr) != size_fin){
         printf("Error reading bv data into temp_noise\n");
            fprintf(fptr,"Error reading bv data into temp_noise\n");
        }


		fclose(fptr);
        
		for (i=0;i<size_fin;i++){
			noise_array[i] += temp_noise[i];
		}


		free(temp_noise);

		// write_out data

		if ((fptr = fopen(filename_real,"w")) != NULL){
			if (fwrite(noise_array,sizeof(double complex),size_fin,fptr) != size_fin){
				fprintf(flog,"Failed to increment observation noise file with incorrect number of members: %s\n",filename_real);
				fflush(flog);
				printf("Failed to increment observation noise file with incorrect number of members: %s\n",filename_real);
			} else {
			fclose(fptr);
			}
		
		} else {
		fprintf(flog,"Failed to increment observation noise file: %s\n",filename_real);
		}


	}


	if ((fptrv = fopen(filename_real2,"r")) == NULL){

		if ((fptrv = fopen(filename_real2,"w")) != NULL){
			fwrite(bdag_v_real,sizeof(double complex),size_fin,fptrv);
			fclose(fptrv);
		} else {
		fprintf(flog,"Failed to write observation visibility file: %s\n",filename_real2);
		printf("Failed to write observation visibility file: %s\n",filename_real2);
		}

	} else {
	

		if (temp_v_real == NULL) temp_v_real = malloc(size_fin*sizeof(double complex));
		if (fread(temp_v_real,sizeof(double complex),size_fin,fptrv) != size_fin) printf("Error reading bv data into temp_v_real\n");

		fclose(fptrv);

		for (i=0;i<size_fin;i++){
			bdag_v_real[i] += temp_v_real[i];
		}

		free(temp_v_real);

		// write_out data

		if ((fptrv = fopen(filename_real2,"w")) != NULL){
			if (fwrite(bdag_v_real,sizeof(double complex),size_fin,fptrv) != size_fin){
				fprintf(flog,"Failed to increment observation vis file with incorrect number of members: %s\n",filename_real2);
				fflush(flog);
				printf("Failed to increment observation vis file with incorrect number of members: %s\n",filename_real2);
			} else {
			fclose(fptrv);
			}
		} else {
		fprintf(flog,"Failed to increment observation visibility file: %s\n",filename_real2);
		}


	}	


        
        
        /* DIFF FILES */
        
        // Check to see if the file already exists
        
        if ((fptr = fopen(filename_reald,"r")) == NULL){
            
            if ((fptr = fopen(filename_reald,"w")) != NULL){
                fwrite(noisediff_array,sizeof(double complex),size_fin,fptr);
                fclose(fptr);
            } else {
                fprintf(flog,"Failed to write observation noise file: %s\n",filename_reald);
            }
            
        } else {
            
            //		printf("File exists... opening...\n");
            
            temp_noisediff = malloc(size_fin*sizeof(double complex));
            assert(temp_noisediff != NULL);
            if (fread(temp_noisediff,sizeof(double complex),size_fin,fptr) != size_fin){
                printf("Error reading bv data into temp_noise\n");
                fprintf(fptr,"Error reading bv data into temp_noise\n");
            }
            
            
            fclose(fptr);
            
            for (i=0;i<size_fin;i++){
                noisediff_array[i] += temp_noisediff[i];
            }
            
            
            free(temp_noisediff);
            
            // write_out data
            
            if ((fptr = fopen(filename_reald,"w")) != NULL){
                if (fwrite(noisediff_array,sizeof(double complex),size_fin,fptr) != size_fin){
                    fprintf(flog,"Failed to increment observation noise file with incorrect number of members: %s\n",filename_reald);
                    fflush(flog);
                    printf("Failed to increment observation noise file with incorrect number of members: %s\n",filename_reald);
                } else {
                    fclose(fptr);
                }
                
            } else {
                fprintf(flog,"Failed to increment observation noise file: %s\n",filename_reald);
            }
            
            
        }
        
        
        if ((fptrv = fopen(filename_real2d,"r")) == NULL){
            
            if ((fptrv = fopen(filename_real2d,"w")) != NULL){
                fwrite(bdag_v_realdiff,sizeof(double complex),size_fin,fptrv);
                fclose(fptrv);
            } else {
                fprintf(flog,"Failed to write observation visibility file: %s\n",filename_real2d);
                printf("Failed to write observation visibility file: %s\n",filename_real2d);
            }
            
        } else {
            
            
            if (temp_v_realdiff == NULL) temp_v_realdiff = malloc(size_fin*sizeof(double complex));
            if (fread(temp_v_realdiff,sizeof(double complex),size_fin,fptrv) != size_fin) printf("Error reading bv data into temp_v_real\n");
            
            fclose(fptrv);
            
            for (i=0;i<size_fin;i++){
                bdag_v_realdiff[i] += temp_v_realdiff[i];
            }
            
            free(temp_v_realdiff);
            
            // write_out data
            
            if ((fptrv = fopen(filename_real2d,"w")) != NULL){
                if (fwrite(bdag_v_realdiff,sizeof(double complex),size_fin,fptrv) != size_fin){
                    fprintf(flog,"Failed to increment observation vis file with incorrect number of members: %s\n",filename_real2d);
                    fflush(flog);
                    printf("Failed to increment observation vis file with incorrect number of members: %s\n",filename_real2d);
                } else {
                    fclose(fptrv);
                }
            } else {
                fprintf(flog,"Failed to increment observation visibility file: %s\n",filename_real2d);
            }
            
            
        }	
        

        
        
        /* WEIGHTS FILE */
        
        /* write-out flags file of sampled uv locations */

	if ((fptrflags = fopen(filename_flags,"r")) == NULL){

		if ((fptrflags = fopen(filename_flags,"w")) != NULL){
			fwrite(flag_array,sizeof(double complex),size_fin,fptrflags);
			fclose(fptrflags);
		} else {
		fprintf(flog,"Failed to write observation flags file: %s\n",filename_flags);
		}

	} else {

		temp_flags = malloc(size_fin*sizeof(double complex));
		assert(temp_flags != NULL);
		if (fread(temp_flags,sizeof(double complex),size_fin,fptrflags) != size_fin) printf("Error reading flags data into temp_flags\n");

		fclose(fptrflags);

		for (i=0;i<size_fin;i++){
			flag_array[i] += temp_flags[i];
		}

		free(temp_flags);

		// write_out data

		if ((fptrflags = fopen(filename_flags,"w")) != NULL){
			fwrite(flag_array,sizeof(double complex),size_fin,fptrflags);
			fclose(fptrflags);
		} else {
		fprintf(flog,"Failed to increment observation flags file: %s\n",filename_flags);
		}

	}



	}

/*************** FREE ARRAYS *********************/

	free(noise_array);
	free(bdag_v_real);
	free(flag_array);
    free(noisediff_array);
    free(bdag_v_realdiff);


return 0;

}





/******************************
 ! NAME:		init_pb
 ! PURPOSE:		set-up PB to be used - uses pre-defined files when ALT/AZ used
 ! ARGUMENTS:	
 ! RETURNS:		integer (0=success)
******************************/

beamdataf *init_pbf(int freq_index, int pol) {
 beamdataf *obj=NULL;
 FILE *fp=NULL,*fpu=NULL;
 int ii,j,k,l,ch,r;
 char beamfilename[1024],ufile[1024];
 float norm=0.;
//    float **normal=NULL,**total_norm_sq=NULL;

  /* allocate space for data structure */

  obj= calloc(1,sizeof(beamdataf));
  if(obj==NULL) {
    fprintf(stderr,"primary_beam: no malloc for main struct\n");
    return NULL;
  }

   /* read-in Fourier-transformed beam file for this frequency channel and polarization */
	obj->real = create6DMatrixFloat(2,NUM_FINE,NUM_W_PLANES,NUM_POINT,SIZE_BEAM,SIZE_BEAM);
	obj->imag = create6DMatrixFloat(2,NUM_FINE,NUM_W_PLANES,NUM_POINT,SIZE_BEAM,SIZE_BEAM);
	/* No w-splitting settings
  obj->real = create3DMatrixFloat(NUM_POINT,SIZE_BEAM,SIZE_BEAM);
  obj->imag = create3DMatrixFloat(NUM_POINT,SIZE_BEAM,SIZE_BEAM);
	*/

  obj->u_tab = calloc(SIZE_BEAM,sizeof(float));
    obj->normal = create3DMatrixFloat(NUM_FINE,NUM_W_PLANES,NUM_POINT);
    obj->total_norm_sq = create3DMatrixFloat(NUM_FINE,NUM_W_PLANES,NUM_POINT);


    /* original RBW beam */
    
    //     if (pol == 0) sprintf(beamfilename,"%sbeam_ft_%03d_9_15_xx_sq.dat",getenv("BEAMDIR"),freq_index);
    //     if (pol == 1) sprintf(beamfilename,"%sbeam_ft_%03d_9_15_yy_sq.dat",getenv("BEAMDIR"),freq_index);
    
    /* New Curtin beam with fine channel linear interpolation!! */
    
		    if (pol == 0) sprintf(beamfilename,"%sbeam_ft_%03d_9_15_xx_finechan.dat",getenv("BEAMDIR"),freq_index);
		    if (pol == 1) sprintf(beamfilename,"%sbeam_ft_%03d_9_15_yy_finechan.dat",getenv("BEAMDIR"),freq_index);
   
	//	sprintf(beamfilename,"%sbeam_ft_%03d_9_15_all_gauss_fine.dat",getenv("BEAMDIR"),freq_index);

    
 //   if (pol == 0) sprintf(beamfilename,"%sbeam_ft_%03d_9_15_xx_sq_exp3.dat",getenv("BEAMDIR"),freq_index);
 //   if (pol == 1) sprintf(beamfilename,"%sbeam_ft_%03d_9_15_yy_sq_exp3.dat",getenv("BEAMDIR"),freq_index);
    
		/* No w-splitting settings
		if (pol == 0) sprintf(beamfilename,"%sbeam_ft_%03d_pt%01d_xx.dat",getenv("BEAMDIR"),freq_index,NUM_POINT);
		if (pol == 1) sprintf(beamfilename,"%sbeam_ft_%03d_pt%01d_yy.dat",getenv("BEAMDIR"),freq_index,NUM_POINT);
		*/

		printf("filename: %s\n",beamfilename);

		if ((fp = fopen(beamfilename, "r")) == NULL){
		fprintf(stderr,"primary_beam: cannot open beam file\n");
		return NULL;
		}
		for (r=0;r<2;r++){
        for (ch=0;ch<NUM_FINE;ch++){
		for (k=0;k<NUM_W_PLANES;k++){	
		for (l=0;l<NUM_POINT;l++){


  		for (ii=0;ii<SIZE_BEAM;ii++){
    		for (j=0;j<SIZE_BEAM;j++){
      			if (fread(&obj->real[r][ch][k][l][j][ii],sizeof(float),1,fp) != 1) printf("Error reading beam file.\n");
      			if (fread(&obj->imag[r][ch][k][l][j][ii],sizeof(float),1,fp) != 1) printf("Error reading beam file.\n");

    		}
  		}

        }
		}
		}
		}
		fclose(fp);

    for (ch=0;ch<NUM_FINE;ch++){
        for (k=0;k<NUM_W_PLANES;k++){
            for (l=0;l<NUM_POINT;l++){
    
                for (ii=0;ii<SIZE_BEAM;ii++){
                        for (j=0;j<SIZE_BEAM;j++){
        
                obj->normal[ch][k][l] += sqrt((obj->real[0][ch][k][l][j][ii])*(obj->real[0][ch][k][l][j][ii]) + (obj->imag[0][ch][k][l][j][ii])*(obj->imag[0][ch][k][l][j][ii]));
               // /(DELTA_U/INTRINSIC_DELTA_U)/(DELTA_U/INTRINSIC_DELTA_U);
                obj->total_norm_sq[ch][k][l] += sqrt((obj->real[1][ch][k][l][j][ii])*(obj->real[1][ch][k][l][j][ii]) + (obj->imag[1][ch][k][l][j][ii])*(obj->imag[1][ch][k][l][j][ii]));
               //  /(DELTA_U/INTRINSIC_DELTA_U)/(DELTA_U/INTRINSIC_DELTA_U);
            
                        }
                }
                                                         
            }
        }
    }
    
    
		/* Introduce new normalisation of beams - March 16, 2015 */

		sprintf(ufile,"%sbeam_u_%d.dat",getenv("BEAMDIR"),SIZE_BEAM);
//		sprintf(ufile,"%sbeam_u_33.dat",getenv("BEAMDIR"));
//        sprintf(ufile,"%sbeam_u_151.dat",getenv("BEAMDIR"));
		if ((fpu = fopen(ufile , "r")) == NULL){
		fprintf(stderr,"primary_beam: cannot open beam u file\n");
		return NULL;
		}
		for (j=0;j<SIZE_BEAM;j++) if(fread(&obj->u_tab[j],sizeof(float),1,fpu) != 1) printf("Error reading beam spacing file.\n");
		fclose(fpu);


return obj;

}


/******************************
 ! NAME:		init_pb_delays
 ! PURPOSE:		set-up PB to be used - uses beamformer delays to compute beam pointing
 ! ARGUMENTS:	pointer to delay array
 ! RETURNS:		integer (0=success)
******************************/

int init_pb_delays(float delays) {




return 0;


}


double modd(double val){

    if ((val > 0.)&&(val < M_PI)) return val;
    if ((val < 0.)&&(val > -M_PI)) return val;
    if ((val > 0.)&&(val >= M_PI)) return val - 2.*M_PI;
    if ((val < 0.)&&(val <= -M_PI)) return val + 2.*M_PI;
    
    

}




/******************************
 ! NAME:		get_pbsq_values_fine
 ! PURPOSE:		get primary beam values at the specified uv points in the grid - fine frequency channels
 ! ARGUMENTS:		pointer to delay array
 ! RETURNS:		integer (0=success)
 ******************************/

int get_pbsq_values_fine(beamdataf *bdata, int w_index, int point, int kernel_size, double *uvalues, double *vvalues, double uu, double vv, double *beamsq_real, double *beamsq_imag, float factor_scale, int ch) {
    
    beam_diff_mwa_vector_2D(bdata->real[0][ch][w_index][point],bdata->imag[0][ch][w_index][point],bdata->u_tab,uvalues,uu,vvalues,vv,kernel_size,beamsq_real,beamsq_imag,factor_scale);
    
    /* No w-splitting */
    //beam_diff_mwa_vector_2D(bdata->real[point],bdata->imag[point],bdata->u_tab,uvalues,uu,vvalues,vv,kernel_size,beam_real,beam_imag);
    
    return 0;
    
    
}




/******************************
 ! NAME:		free_pb
 ! PURPOSE:		frees the beam data structure
 ! ARGUMENTS:		beamdata structure type
 ! RETURNS:		void
******************************/


void free_pbf(beamdataf *b){

 int j,k,l,r,ch;

if (b->u_tab !=NULL) free(b->u_tab);

for (r=0;r<2;r++){
for (ch=0;ch<NUM_FINE;ch++){
for (k=0;k<NUM_W_PLANES;k++){
	for (l=0;l<NUM_POINT;l++){
		for (j=0;j<SIZE_BEAM;j++){
			free(b->real[r][ch][k][l][j]);
			free(b->imag[r][ch][k][l][j]);
	//		free(b->real[l][j]);
	//		free(b->imag[l][j]);

        }
		}

	}
}
}
    
    for (k=0;k<NUM_W_PLANES;k++){
        for (l=0;l<NUM_POINT;l++){
    free(b->normal[k][l]);
    free(b->total_norm_sq[k][l]);
        }
        }
    
free(b);


}



/******************************
 ! NAME:		beam_diff_mwa_vector_2D
 ! PURPOSE:		computes the actual beam values at the interpolated points
 ! ARGUMENTS:	high-res beam real and imaginary, returns interpolated real and imaginary beam values
 ! RETURNS:		int (0=success)
******************************/

int beam_diff_mwa_vector_2D(float **beam_r, float **beam_i, float *u_tab, double *u,double uprime,double *v,double vprime,int size,double *out_real,double *out_imag, float factor_scale){

  double xloc,yloc;
  int i;
   float utemp=0.,vtemp=0.;
   float sqrt_size=1.;

sqrt_size = sqrt((float) size);

  /* interpolate beam at given uprime,vprime values */
//    factor_scale=1.;
for (i=0;i<size;i++){
    utemp = u[i];///factor_scale;
    vtemp = v[i];///factor_scale;
	if ((sqrt((utemp-uprime)*(utemp-uprime)) >= u_tab[0])||(sqrt((vtemp-vprime)*(vtemp-vprime)))) {
	  xloc = ((utemp-uprime)*factor_scale-u_tab[0])/(u_tab[1]-u_tab[0]);
	  yloc = ((vtemp-vprime)*factor_scale-u_tab[0])/(u_tab[1]-u_tab[0]);


//if ((yloc < 4)||(xloc < 4)||(yloc > sqrt_size-4)||(xloc > sqrt_size-4)){
	interp_bilin(beam_r,xloc,yloc,&out_real[i]);
	interp_bilin(beam_i,xloc,yloc,&out_imag[i]);

        /*
} else {

	interp(beam_r,xloc,yloc,&out_real[i]);
	interp(beam_i,xloc,yloc,&out_imag[i]);

}
         */


//	printf("i %d u %lg uprime %lg utab0 %lg xloc %lg v %lg vprime %lg yloc %lg interp %lg %lg %lg\n",i,utemp,uprime,u_tab[0],xloc,vtemp,vprime,yloc,out_imag[i],beam_r[0][0],beam_r[1][0]);
//	printf("i %d u %lg uprime %lg xloc %lg v %lg vprime %lg yloc %lg interp %lg\n",i,u[i],uprime,xloc,v[i],vprime,yloc,out_imag[i]);
	} else {
	if (debug) printf("off edge u[i] %g uprime %g\n",u[i],uprime);
	out_real[i] = 0.;
	out_imag[i] = 0.;
	}
}


return 0;
}


/***************************************************
 ! NAME:		interp_bilin
 ! PURPOSE:		Bilinear interpolation for complex numbers
 sx,sy: size of grid x and y axes
 u,v: desired location to interpolate to in pixel units.
 
 ***************************************************/

int interp_bilin(float **uvgrid, double u, double v, double *res_re) {
    int m,n;
    float p,q;
    float Q11,Q12,Q21,Q22;
    float re=0;
    
    m = (int)u;   /* nearest integer cell number such that offset is positive */
    n = (int)v;
    if (u < 0) m -= 1;
    if (v < 0) n -= 1;
    p = u-m;	 /* p,q become positive offsets within the cell */
    q = v-n;
    
    Q11 = uvgrid[m][n];
    Q21 = uvgrid[m+1][n];
    Q12 = uvgrid[m][n+1];
    Q22 = uvgrid[m+1][n+1];
    
    re = Q11*(1.-p)*(1.-q) + Q21*p*(1.-q) + Q12*(1.-p)*q + Q22*p*q;
    
 //   printf("bilinear xloc %f yloc %f uvgrid %f %f %f %f output %f\n",u,v,uvgrid[m][n],uvgrid[m+1][n],uvgrid[m][n+1],uvgrid[m+1][n+1],re);
    *res_re = re;
    
    return 0;
}


/************************************************/

/***************************************************
 ! NAME:		interp
 ! PURPOSE:		2D cublic spline interpolation for complex numbers
   sx,sy: size of grid x and y axes
   u,v: desired location to interpolate to in pixel units.

***************************************************/

#define SIZ_SPLINE 7

int interp(float **uvgrid, double u, double v, double *res_re) {
  int m,n,i,j;
  float p,q;
  float x1[SIZ_SPLINE],x2[SIZ_SPLINE];
  static float *y[SIZ_SPLINE],*y2[SIZ_SPLINE];
  static int init=0;
  float re=0;

  m = (int)u;   /* nearest integer cell number such that offset is positive */
  n = (int)v;
  if (u < 0) m -= 1;
  if (v < 0) n -= 1;
  p = u-m;	 /* p,q become positive offsets within the cell */
  q = v-n;

 //   if (debug) {
//   printf("interp_vis: p,q,m,n: %g,%g,%d,%d \n",p,q,m,n);
 /*   for(j=0; j<SIZ_SPLINE; j++) {
      for(i=0; i<SIZ_SPLINE; i++) {
	val = getUVCell2D(i+m-1,j+n-1,sx,sy, uvgrid);
        printf("(%g,%g) ",creal(val),cimag(val));
      }
      printf("\n");
    } */
//  }
  

  /* initialise the x and y arrays for cubic spline interp */
    for (i=0; i<SIZ_SPLINE; i++){
        x1[i]=i-1.0; /* coords of columns (V) from UV plane */
        x2[i]=i-1.0; /* coords of rows (U) from UV plane */
    }
    
    
  if (!init) {
    for (i=0; i<SIZ_SPLINE; i++){
      y[i]  = calloc(SIZ_SPLINE,sizeof(float));
      y2[i] = calloc(SIZ_SPLINE,sizeof(float));
    }
    init=1;
    
  }

  re=0.0;

    /* set y values for spline (REAL) 
    for (j=0; j<SIZ_SPLINE; j++){          loop over rows 
      for (i=0; i<SIZ_SPLINE; i++){        loop over cols 
	val=getUVCell2D(m+i-1,n+j-1,sx,sy, uvgrid);
        y[j][i] = creal(val);
      }
    }  
*/

//printf("Here\n");

    /* set y values for spline (REAL) */
    for (j=0; j<SIZ_SPLINE; j++){         /* loop over rows */
      for (i=0; i<SIZ_SPLINE; i++){       /* loop over cols */
//printf("%d %d %d %d\n",m,n,j,i);
	y[j][i]=uvgrid[m+i-1][n+j-1];
      }
    }

//printf("Here2a\n");

    /* initialise spline coeffs */
    splie2(x1,x2,y,SIZ_SPLINE,SIZ_SPLINE,y2);
    /* now do the interp (REAL) */
    splin2(x1,x2,y,y2,SIZ_SPLINE,SIZ_SPLINE,q,p,&re);

//printf("Here2\n");
    
//     printf("spline xloc %f yloc %f uvgrid %f %f %f %f output %f\n",u,v,uvgrid[m][n],uvgrid[m+1][n],uvgrid[m][n+1],uvgrid[m+1][n+1],re);

  *res_re = re;
    
//printf("Here3\n");
 //   free(x1);
 //   free(x2);
 //   free(y);
  //  free(y2);
    
//printf("Here4\n");

  return 0;
}

double gauss_noise (/* value1, value2) */
double value1){

double fac, rsq, v1, v2;

while (1)
	{
									// Random location in unit square
	v1 = 2.0 * drand48() - 1.0;
	v2 = 2.0 * drand48() - 1.0;
	rsq = v1*v1 + v2*v2;
                                        // Take only those inside unit circle
	if (rsq < 1.0) break;
	}
                                        // Following the recipe ...
fac = sqrt (-2.0 * log (rsq) / rsq);

return v1*fac;

}


/**************************************
 ! NAME:		create4DMatrix
 ! PURPOSE:		create a 4D double array
 ! RETURNS:		pointer to a 4D Float array
**************************************/

float ****create4DMatrixFloat(int sizex, int sizey, int sizez, int sizew){
float ****output;
 int k,l,m;

 output = calloc(sizex,sizeof(float***));
 for (k=0;k<sizex;k++){
   output[k] = calloc(sizey,sizeof(float**));
 }
 for (k=0;k<sizex;k++){
   for (l=0;l<sizey;l++){
     output[k][l] = calloc(sizez,sizeof(float*));
   }
 }
 for (k=0;k<sizex;k++){
   for (l=0;l<sizey;l++){
     for (m=0;m<sizez;m++){
       output[k][l][m] = calloc(sizew,sizeof(float));
     }
   }
 }
 return output;


}


/***************************************************/


/**************************************
 ! NAME:		create5DMatrix
 ! PURPOSE:		create a 5D double array
 ! RETURNS:		pointer to a 5D Float array
**************************************/

float *****create5DMatrixFloat(int sizer, int sizex, int sizey, int sizez, int sizew){
float *****output;
 int k,l,m,r;

output = calloc(sizer,sizeof(float****));
 for (r=0;r<sizer;r++){
   output[r] = calloc(sizex,sizeof(float***));
 }
 for (r=0;r<sizer;r++){
 for (k=0;k<sizex;k++){
   output[r][k] = calloc(sizey,sizeof(float**));
 }
}
for (r=0;r<sizer;r++){
 for (k=0;k<sizex;k++){
   for (l=0;l<sizey;l++){
     output[r][k][l] = calloc(sizez,sizeof(float*));
   }
 }
}
for (r=0;r<sizer;r++){
 for (k=0;k<sizex;k++){
   for (l=0;l<sizey;l++){
     for (m=0;m<sizez;m++){
       output[r][k][l][m] = calloc(sizew,sizeof(float));
     }
   }
 }
}
 return output;


}


/***************************************************/

/**************************************
 ! NAME:		create6DMatrix
 ! PURPOSE:		create a 6D double array
 ! RETURNS:		pointer to a 6D Float array
 **************************************/

float ******create6DMatrixFloat(int sizech, int sizer, int sizex, int sizey, int sizez, int sizew){
    float ******output;
    int k,l,m,r,ch;
    
    output = calloc(sizech,sizeof(float*****));
    for (ch=0;ch<sizech;ch++){
        output[ch] = calloc(sizer,sizeof(float****));
    }
    
    for (ch=0;ch<sizech;ch++){
        for (r=0;r<sizer;r++){
            output[ch][r] = calloc(sizex,sizeof(float***));
        }
    }
    
    for (ch=0;ch<sizech;ch++){
    for (r=0;r<sizer;r++){
        for (k=0;k<sizex;k++){
            output[ch][r][k] = calloc(sizey,sizeof(float**));
        }
    }
    }
    for (ch=0;ch<sizech;ch++){
    for (r=0;r<sizer;r++){
        for (k=0;k<sizex;k++){
            for (l=0;l<sizey;l++){
                output[ch][r][k][l] = calloc(sizez,sizeof(float*));
            }
        }
    }
        }
    for (ch=0;ch<sizech;ch++){
    for (r=0;r<sizer;r++){
        for (k=0;k<sizex;k++){
            for (l=0;l<sizey;l++){
                for (m=0;m<sizez;m++){
                    output[ch][r][k][l][m] = calloc(sizew,sizeof(float));
                }
            }
        }
    }
    }
    return output;
    
    
}


/***************************************************/





/**************************************
 ! NAME:		create3DMatrix
 ! PURPOSE:		create a 3D double array
 ! RETURNS:		pointer to a 3D Float array
**************************************/

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

/************************************************/



