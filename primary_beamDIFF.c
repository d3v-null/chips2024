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
#include <star/pal.h>
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

int calc_bv_flags_diff(uvdata *data, uvdata *data2, beamdata *beam_data, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext, int *flag14m){

FILE *fptr=NULL,*fptrv=NULL,*fptrflags=NULL,*flog=NULL;
int flag=0,flag2=0,w_index,obs_year,obs_month,obs_day,i,j,k,ii,jj,xx,beam_size,kernel_size,top,*include_vis=NULL;
double w_centre,checksum_real=0.0,checksum_imag=0.0,checksum_abs=0.0,checksum_cent=0.0,v1,v2,v3,v4;
long contrib_vis=0,total_vis=0;
float vis_r1,vis_r2,vis_i1,vis_i2;
char filename_real[1024],filename_real2[1024],filename_reald[1024],filename_real2d[1024],filename_flags[1024],debugfile[1024];
double complex  *bdag_v_real=NULL,*bdag_v_realdiff=NULL,*temp_v_real=NULL,*temp_noise=NULL,*temp_v_realdiff=NULL,*temp_noisediff=NULL,*noise_array=NULL,*noisediff_array=NULL;
double *flag_array=NULL,*temp_flags=NULL,*flagw_array=NULL,*temp_flagsw=NULL;
float vis_rtot,vis_rdiff,vis_itot,vis_idiff,w_temp=0.,kernel_size_float,visnoise,weight,normal=0.,total_norm_sq=0.,frequency_lower=0.,factor_scale=1.,l8,additup=0.;
double distance=0.,*beam_real=NULL,*beam_imag=NULL,*u_lex_small=NULL,*v_lex_small=NULL,uu,vv,ww,norm,*beamsq_real=NULL,*beamsq_imag=NULL,*beamsq_real2=NULL,*beamsq_imag2=NULL;


//    printf("internal flag14m %d %d %d %d\n",flag14m[0],flag14m[10],flag14m[4567],flag14m[8127]);

    
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
norm = 1./(DELTA_U*DELTA_U/intrinsic_res/intrinsic_res);
printf("Beam normalization: %g\n",norm);
 
visnoise = 1.;
//norm2=1.;
    
/* define frequency at low end of coarse channel */
    
    frequency_lower = freq_index*COARSE_CHAN_WIDTH + LOWER_FREQ;
    factor_scale = frequency/frequency_lower;
//    factor_scale = 1.;
    printf("factor_scale %f\n",factor_scale);

/* Define dense vector for bdag_v accumulation and flags for uv-sampling */

	bdag_v_real = calloc(size_fin,sizeof(double complex));
	noise_array = calloc(size_fin,sizeof(double complex));
    bdag_v_realdiff = calloc(size_fin,sizeof(double complex));
    noisediff_array = calloc(size_fin,sizeof(double complex));
	flag_array = calloc(size_fin,sizeof(double));


    
    l8 = global_period/86400.*2.*M_PI;     // extra rotation for second step
 //   l8 = 0.;
    
	// set weights to zero for 14m east-west baselines
    
 //   printf("internal flag14m %d %d %d %d\n",flag14m[0],flag14m[10],flag14m[4567],flag14m[8127]);
    
//	int ws = 10 - wstack;
//	wstack = ws;

for (wstack=0;wstack<NUM_W_STACK;wstack++){
//  for (wstack=0;wstack<1;wstack++){
    
	w_centre = (wstack*DELTA_WSTACK);
    
//	printf("wstack %d w_centre %lg\n",wstack,w_centre);

   /*****************************************************************/
   /* Loop over baselines - first loop to compute memory allocation */

	include_vis = calloc(data->n_baselines[0],sizeof(int));

	for (i=0;i<data->n_baselines[0];i++){
     //    printf("baseline: %d\n",i);
				
     /* Compute whether uv point lies within range of uv-grid and wstack */

		distance = sqrt(data->u[0][i]*data->u[0][i]+data->v[0][i]*data->v[0][i])*(frequency);
	//	if (debug) fprintf(flog,"i %d distance: %lg w: %g, weight %g, WMAX %f, wcentre %f\n",i,distance,data->w[0][i]*frequency,data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)],WMAX,w_centre);
        

 
 		if ((distance >= UMAX-SIZE_BEAM/2*INTRINSIC_DELTA_U) || (distance < 0.5) || (sqrt(data->w[0][i]*frequency*data->w[0][i]*frequency) >= WMAX) || ((sqrt((sqrt(data->w[0][i]*data->w[0][i])*frequency-w_centre)*(sqrt(data->w[0][i]*data->w[0][i])*frequency-w_centre))) > WSTACK_DELT_TOT/2.) || (data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] <= 0.) || (flag14m[i] == 0)){
            

		} else {

			contrib_vis++;
			include_vis[i] = 1;

		}




	}



	if (debug) fprintf(flog,"contrib vis: %ld\n",contrib_vis);

	/* set-up beam to contain enough entries to comfortably contain beam
			  region */
	
//	beam_size = (int) BEAM_SIZE_FLOAT*(1.+w_centre/10.);
    beam_size = (int) BEAM_SIZE_FLOAT;
 //   printf("beam size %d\n",beam_size);
		
	if (((float)beam_size)/2. == (int)((float)(beam_size)/2.)) beam_size = beam_size+1;  /* make odd */

	kernel_size_float = beam_size*beam_size;
	kernel_size = (int)kernel_size_float;

//	int size_beam = (int)((float)beam_size*(float)beam_size*(float)beam_size*(float)beam_size);

   /*****************************************************************/
   /* Loop over baselines - second loop to actually process data */

	if (contrib_vis == 0){
	} else {

	for (i=0;i<data->n_baselines[0];i++){
     //    printf("baseline: %d\n",i);
				
 //       if (data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] > 64.) printf("weight large!\n");

		if (include_vis[i] != 1){
		} else {

			//contrib_vis++;

			uu = data->u[0][i]*(frequency);
			vv = data->v[0][i]*(frequency);
			ww = data->w[0][i]*(frequency);
			flag = 0;

			if (((HALF_PLANE_FLAG)&&(uu < 0.0))||(ww < 0.)){
   // test         if (((HALF_PLANE_FLAG)&&(uu < 0.0))){
			flag = 1;
			uu = -uu;
			vv = -vv;
			ww = -ww;
			} 

			/* remove points in -ve v for u=0 */
			if ((HALF_PLANE_FLAG)&&(round(uu/DELTA_U) == 0)&&(vv < 0.)){
			flag = 1;
			vv = -vv;
			} 

    			if (debug) fprintf(flog,"u %g v %g w %g\n",uu,vv,ww);

			if (debug) fflush(flog);

			u_lex_small = calloc(kernel_size,sizeof(double));
			v_lex_small = calloc(kernel_size,sizeof(double));
			long loc = 0;
//      		printf("beam size %d kernel size %d\n",beam_size,kernel_size);

    			   /* lexicographic u and v co-ordinates of small *sky* Fourier plane around the Fourier location */
			for (k=0;k<beam_size;k++){
				for (j=0;j<beam_size;j++){
								 
					u_lex_small[loc] = (double) ((k)-(int) beam_size/2)*DELTA_U + round(uu/DELTA_U)*DELTA_U;
					v_lex_small[loc] = (double) ((j)-(int) beam_size/2)*DELTA_U + round(vv/DELTA_U)*DELTA_U;
					if (debug) printf("loc %ld u_lex %g v_lex %g\n",loc,u_lex_small[loc],v_lex_small[loc]);
					loc++;
				}
			}
								   
			beam_real = calloc(kernel_size,sizeof(double));
			beam_imag = calloc(kernel_size,sizeof(double));
   
			beamsq_real = calloc(kernel_size,sizeof(double));
			beamsq_imag = calloc(kernel_size,sizeof(double));
            
            beamsq_real2 = calloc(kernel_size,sizeof(double));
            beamsq_imag2 = calloc(kernel_size,sizeof(double));

			/* GET PRIMARY BEAM VALUES HERE */

            get_pbsq_values(beam_data,wstack,point,kernel_size,u_lex_small,v_lex_small,uu,vv,beamsq_real,beamsq_imag,1.);
//			get_pbsq_values(beam_data,wstack,point,kernel_size,u_lex_small,v_lex_small,uu,vv,beam_real,beam_imag,1.);

            total_norm_sq = beam_data->total_norm_sq[wstack][point];
//		total_norm_sq = 1.;
            normal = (beam_data->total_norm_sq[wstack][point]);   ///sqrt(factor_scale);
            
            
			// Rotate beam to phase centre
			for (k=0;k<kernel_size;k++){

                float rotat1 = -(2.*M_PI*(u_lex_small[k]-uu)*sin(modd(l_rot)));
                float rotat2 = -(2.*M_PI*(v_lex_small[k]-vv)*sin(modd(m_rot)));
                
                float rotat3 = -(2.*M_PI*(u_lex_small[k]-uu)*sin(modd(l_rot-l8)));

                
                //				printf("Rotation to beam centre (rad): %g\n",rotat);
                float temp1 = beamsq_real[k]*cos(rotat1) + beamsq_imag[k]*sin(rotat1);
                float temp2 = beamsq_imag[k]*cos(rotat1) - beamsq_real[k]*sin(rotat1);

                
                float temp3 = beamsq_real[k]*cos(rotat3) + beamsq_imag[k]*sin(rotat3);
                float temp4 = beamsq_imag[k]*cos(rotat3) - beamsq_real[k]*sin(rotat3);

                
                // testing rotation direction
                //				float temp1 = beam_real[k]*cos(rotat) - beam_imag[k]*sin(rotat);
                //				float temp2 = beam_imag[k]*cos(rotat) + beam_real[k]*sin(rotat);
                
                beamsq_real[k] = temp1;
                beamsq_imag[k] = temp2;
                
                beamsq_real2[k] = temp3;
                beamsq_imag2[k] = temp4;
                

                
                
                float temp11 = beamsq_real[k]*cos(rotat2) + beamsq_imag[k]*sin(rotat2);
                float temp22 = beamsq_imag[k]*cos(rotat2) - beamsq_real[k]*sin(rotat2);
 
                float temp33 = beamsq_real2[k]*cos(rotat2) + beamsq_imag2[k]*sin(rotat2);
                float temp44 = beamsq_imag2[k]*cos(rotat2) - beamsq_real2[k]*sin(rotat2);

                
                // testing rotation direction
                //				float temp1 = beam_real[k]*cos(rotat) - beam_imag[k]*sin(rotat);
                //				float temp2 = beam_imag[k]*cos(rotat) + beam_real[k]*sin(rotat);
                
                beamsq_real[k] = temp11;
                beamsq_imag[k] = temp22;
                

                beamsq_real2[k] = temp33;
                beamsq_imag2[k] = temp44;


      //   			  printf("k %d beamr %g beam i %g\n",k,beam_real[k],beam_imag[k]);
			}
	
			/* visibility data from uvfits input file -- cc if in -ve u plane
			factor of pol is to jump to the YY pol after the XX */
/*			vis_r = data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)];
			if (flag == 0){ vis_i = data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1]; } else { vis_i = -data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1];}
 
*/
            
            
        /* ONLY GRID THESE POINTS IF THE WEIGHTS MATCH */
            
 //       if ((data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] == data2->weightdata[0][(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)])&&(data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] == 8.)) {
            
            if ((data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] == data2->weightdata[0][(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)])&&(data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] > 0.)){
           
 /*
			vis_rtot = (data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] + data2->visdata[0][2*(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)]);
			if (flag == 0){ vis_itot = (data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1] + data2->visdata[0][2*(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)+1]); } else { vis_itot = (-data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1]-data2->visdata[0][2*(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)+1]);}

            vis_rdiff = (data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] - data2->visdata[0][2*(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)]);
            if (flag == 0){ vis_idiff = (data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1] - data2->visdata[0][2*(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)+1]); } else { vis_idiff = (-data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1] + data2->visdata[0][2*(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)+1]);}
 */
            
            
            vis_r1 = data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)];
            vis_r2 = data2->visdata[0][2*(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)];
            vis_i1 = data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)+1];
            vis_i2 = data2->visdata[0][2*(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)+1];
  
 
/************************* TESTING ***********************/

//vis_r=1.;
//vis_i=0.;


            v1 = gauss_noise(1.)/sqrt(2.);
            v2 = gauss_noise(1.)/sqrt(2.);
            v3 = gauss_noise(1.)/sqrt(2.);
            v4 = gauss_noise(1.)/sqrt(2.);
  //          printf("v1 %f v2 %f\n",v1,v2);

			/* weights for each visibility */
			weight = 0.5*(data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] + data2->weightdata[0][(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)]);

//printf("ch %d tot %lf %lf diff %lf %lf weight %lf\n",ch,vis_r1,vis_i1,vis_r2,vis_i2,weight);

		if ((debug)&&(weight != 0.0)) fprintf(flog,"chan_width %f period %f weight %f\n",global_chanwidth,global_period,weight);


            weight = sqrt(weight/(global_chanwidth/1.e4*global_period/8.));
	//                weight = 1.;
//                printf("Weight: %f\n",weight);

            
//            printf("pol %d abs vis %f\n",pol,vis_r*vis_r+vis_i*vis_i);

            if ((debug)&&(weight != 0.0)) fprintf(flog,"chan_width %f period %f ch %d pol %d weight %f vis_r %f beam %f loc %d\n",global_chanwidth,global_period,ch,pol,weight,creal(vis_itot),(beamsq_real[10]),xx);

			if (debug) fflush(flog);

		       /* update bdag_b matrix and bdag_v vector according to the footprint of this baseline */

                additup=0.;
                
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

						
					/* MEASURED VISIBILITY DATA CALC */

//					if (debug){
					if ((isnan(vis_rtot) != 0)) printf("NAN NAN NAN NAN %f %f\n",vis_rtot,vis_itot);
//					}

                    /*
					if (flag2 == 0){
					bdag_v_real[xx] += (beam_real[j]*vis_r + beam_imag[j]*vis_i)/normal/norm*weight*weight + I*(beam_real[j]*vis_i - beam_imag[j]*vis_r)/normal/norm*weight*weight;
					} else {
					bdag_v_real[xx] += (beam_real[j]*vis_r + beam_imag[j]*vis_i)/normal/norm*weight*weight - I*(beam_real[j]*vis_i - beam_imag[j]*vis_r)/normal/norm*weight*weight;
					}
*/
                    
                    /*
                    
                    // TOT
                    if (flag2 == 0){
                        bdag_v_real[xx] += (beamsq_real[j]*vis_r1 + beamsq_real2[j]*vis_r2 + beamsq_imag[j]*(-vis_i1)+beamsq_imag2[j]*(-vis_i2))/(total_norm_sq)*weight*weight + I*(beamsq_real[j]*(-vis_i1)+beamsq_real2[j]*(-vis_i2) - (beamsq_imag[j]*vis_r1 + beamsq_imag2[j]*vis_r2))/(total_norm_sq)*weight*weight;
                    } else {
                        bdag_v_real[xx] += (beamsq_real[j]*vis_r1 + beamsq_real2[j]*vis_r2 + beamsq_imag[j]*(-vis_i1)+beamsq_imag2[j]*(-vis_i2))/(total_norm_sq)*weight*weight - I*(beamsq_real[j]*(-vis_i1)+beamsq_real2[j]*(-vis_i2) - (beamsq_imag[j]*vis_r1 + beamsq_imag2[j]*vis_r2))/(total_norm_sq)*weight*weight;
                    }
                    
                    
                    //  DIFF
                    if (flag2 == 0){
                        bdag_v_realdiff[xx] += (beamsq_real[j]*vis_r1 - beamsq_real2[j]*vis_r2 + beamsq_imag[j]*(-vis_i1)-beamsq_imag2[j]*(vis_i2))/(total_norm_sq)*weight*weight + I*(beamsq_real[j]*(-vis_i1)+beamsq_real2[j]*(vis_i2) - (beamsq_imag[j]*vis_r1 - beamsq_imag2[j]*vis_r2))/(total_norm_sq)*weight*weight;
                    } else {
                        bdag_v_realdiff[xx] += (beamsq_real[j]*vis_r1 - beamsq_real2[j]*vis_r2 + beamsq_imag[j]*(-vis_i1)-beamsq_imag2[j]*(vis_i2))/(total_norm_sq)*weight*weight - I*(beamsq_real[j]*(-vis_i1)+beamsq_real2[j]*(vis_i2) - (beamsq_imag[j]*vis_r1 - beamsq_imag2[j]*vis_r2))/(total_norm_sq)*weight*weight;
                    }
                    

                    if (flag == 0){
                        
                        
                        if (flag2 == 0){
                            bdag_v_real[xx] += (beamsq_real[j]*vis_r1 + beamsq_real2[j]*vis_r2 + beamsq_imag[j]*(vis_i1)+beamsq_imag2[j]*(vis_i2))/(total_norm_sq)*weight*weight + I*(beamsq_real[j]*(vis_i1)+beamsq_real2[j]*(vis_i2) - (beamsq_imag[j]*vis_r1 + beamsq_imag2[j]*vis_r2))/(total_norm_sq)*weight*weight;
                        } else {
                            bdag_v_real[xx] += (beamsq_real[j]*vis_r1 + beamsq_real2[j]*vis_r2 + beamsq_imag[j]*(vis_i1)+beamsq_imag2[j]*(vis_i2))/(total_norm_sq)*weight*weight - I*(beamsq_real[j]*(vis_i1)+beamsq_real2[j]*(vis_i2) - (beamsq_imag[j]*vis_r1 + beamsq_imag2[j]*vis_r2))/(total_norm_sq)*weight*weight;
                        }
                        
                        
                        
                        if (flag2 == 0){
                            bdag_v_realdiff[xx] += (beamsq_real[j]*vis_r1 - beamsq_real2[j]*vis_r2 + beamsq_imag[j]*(vis_i1)-beamsq_imag2[j]*(-vis_i2))/(total_norm_sq)*weight*weight + I*(beamsq_real[j]*(vis_i1)+beamsq_real2[j]*(-vis_i2) - (beamsq_imag[j]*vis_r1 - beamsq_imag2[j]*vis_r2))/(total_norm_sq)*weight*weight;
                        } else {
                            bdag_v_realdiff[xx] += (beamsq_real[j]*vis_r1 - beamsq_real2[j]*vis_r2 + beamsq_imag[j]*(vis_i1)-beamsq_imag2[j]*(-vis_i2))/(total_norm_sq)*weight*weight - I*(beamsq_real[j]*(vis_i1)+beamsq_real2[j]*(-vis_i2) - (beamsq_imag[j]*vis_r1 - beamsq_imag2[j]*vis_r2))/(total_norm_sq)*weight*weight;
                        }

                        
                        
                        
                    }
                    */
                    
                    if ((flag == 0)&&(flag2 == 0)){
                        
                        
                        bdag_v_real[xx] += (vis_r1*beamsq_real2[j] + vis_r2*beamsq_real[j] - vis_i1*beamsq_imag2[j] - vis_i2*beamsq_imag[j] )/normal/norm*weight*weight +I*( vis_r1*beamsq_imag2[j] + vis_i1*beamsq_real2[j] + vis_r2*beamsq_imag[j] + vis_i2*beamsq_real[j] )/normal/norm*weight*weight;
                        
                        bdag_v_realdiff[xx] += ( vis_r1*beamsq_real2[j] + vis_i2*beamsq_imag[j] - vis_i1*beamsq_imag2[j] - vis_r2*beamsq_real[j] )/normal/norm*weight*weight +I*( vis_r1*beamsq_imag2[j] + vis_i1*beamsq_real2[j] - vis_r2*beamsq_imag[j] - vis_i2*beamsq_real[j] )/normal/norm*weight*weight;
                        
                        
                    } else {
                        
                        bdag_v_real[xx] += (vis_r1*beamsq_real2[j] + vis_r2*beamsq_real[j] + vis_i1*beamsq_imag2[j] + vis_i2*beamsq_imag[j] )/normal/norm*weight*weight +I*( vis_r1*beamsq_imag2[j] - vis_i1*beamsq_real2[j] + vis_r2*beamsq_imag[j] - vis_i2*beamsq_real[j] )/normal/norm*weight*weight;
                        
                        bdag_v_realdiff[xx] += ( vis_r1*beamsq_real2[j] - vis_i2*beamsq_imag[j] + vis_i1*beamsq_imag2[j] - vis_r2*beamsq_real[j] )/normal/norm*weight*weight +I*( vis_r1*beamsq_imag2[j] - vis_i1*beamsq_real2[j] - vis_r2*beamsq_imag[j] + vis_i2*beamsq_real[j] )/normal/norm*weight*weight;
                        
                        
                    }
                    
                    
                    
                    
                    if ((debug)&&(ch == 9)) fprintf(flog,"bdag %f %f weight %f\n",creal(bdag_v_real[xx]),cimag(bdag_v_real[xx]),weight);
                    
                   
/* EXPECTED NOISE CALC: Compute the expected contribution from a thermal-noise only visibility, located at the same point in the uvw-plane */

                    
                    
			if (flag2 == 0){
			noise_array[xx] += (beamsq_real[j]*(v1+v3) + beamsq_imag[j]*(v2+v4))/normal/norm*weight*weight + I*(-beamsq_imag[j]*(v1+v3) + beamsq_real[j]*(v2+v4))/normal/norm*weight*weight;
			} else {
			noise_array[xx] += (beamsq_real[j]*(v1+v3) + beamsq_imag[j]*(v2+v4))/normal/norm*weight*weight - I*(-beamsq_imag[j]*(v1+v3) + beamsq_real[j]*(v2+v4))/normal/norm*weight*weight;
			}
                    
                    if (flag2 == 0){
                        noisediff_array[xx] += (beamsq_real[j]*(v1-v3) + beamsq_imag[j]*(v2-v4))/normal/norm*weight*weight + I*(-beamsq_imag[j]*(v1-v3) + beamsq_real[j]*(v2-v4))/normal/norm*weight*weight;
                    } else {
                        noisediff_array[xx] += (beamsq_real[j]*(v1-v3) + beamsq_imag[j]*(v2-v4))/normal/norm*weight*weight - I*(-beamsq_imag[j]*(v1-v3) + beamsq_real[j]*(v2-v4))/normal/norm*weight*weight;
                    }
                    

             //       printf("flagarray %lg\n",sqrt(0.5)*sqrt(beamsq_real[j]*beamsq_real[j] + beamsq_imag[j]*beamsq_imag[j] + beamsq_real2[j]*beamsq_real2[j] + beamsq_imag2[j]*beamsq_imag2[j])/total_norm_sq/norm*weight*weight);
                    
		flag_array[xx] += sqrt(0.5)*sqrt(beamsq_real[j]*beamsq_real[j] + beamsq_imag[j]*beamsq_imag[j] + beamsq_real2[j]*beamsq_real2[j] + beamsq_imag2[j]*beamsq_imag2[j])/total_norm_sq/norm*weight*weight;
                     
             //       additup += flag_array[xx];
                    
				}
			}   /* end loop over beam */
                
                
          //      printf("additup %lg %g\n",additup,total_norm_sq);
                
	}    /* end loop over weights being equal */
       

			total_vis++;
			if (debug) fprintf(flog,"total vis %ld\n",total_vis);

            
        
        
			/* clean-up */
			free(beam_real);
			free(beam_imag);
			free(beamsq_real);
			free(beamsq_imag);
            free(beamsq_real2);
            free(beamsq_imag2);
			free(u_lex_small);
			free(v_lex_small);


		}   /* end IF/ELSE loop for baseline inclusion */

	} /* end loop over baselines for a given w stack and pol */
        

	} /* end if statement - don't process anything if no vis contribute!! */

	free(include_vis);

//	printf("Number of baselines contributing %ld of %d, at w-centre %g, polarization %d and frequency %g\n",contrib_vis,data->n_baselines[0],w_centre,pol,frequency);

    } /* End loop over wstacks before writing out files */
    
  //  free(flag14m);

    
	if (debug) fclose(flog);
//	if (debug) exit(1);

	if (contrib_vis > 0){


/************** OUTPUT FILES ********************/
	/* define filenames */

	JD_to_Cal(data->date[0],&obs_year, &obs_month, &obs_day);
	printf("Year %d Month %d Day %d\n",obs_year,obs_month,obs_day);

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
			fwrite(flag_array,sizeof(double),size_fin,fptrflags);
			fclose(fptrflags);
		} else {
		fprintf(flog,"Failed to write observation flags file: %s\n",filename_flags);
		}

	} else {

		temp_flags = malloc(size_fin*sizeof(double));
		assert(temp_flags != NULL);
		if (fread(temp_flags,sizeof(double),size_fin,fptrflags) != size_fin) printf("Error reading flags data into temp_flags\n");

		fclose(fptrflags);

		for (i=0;i<size_fin;i++){
			flag_array[i] += temp_flags[i];
		}

		free(temp_flags);

		// write_out data

		if ((fptrflags = fopen(filename_flags,"w")) != NULL){
			fwrite(flag_array,sizeof(double),size_fin,fptrflags);
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

beamdata *init_pb(int freq_index, int pol) {
 beamdata *obj=NULL;
 FILE *fp=NULL,*fpu=NULL;
 int ii,j,k,l,r;
 char beamfilename[1024],ufile[1024];
 float norm=0.;
//    float **normal=NULL,**total_norm_sq=NULL;

  /* allocate space for data structure */

  obj= calloc(1,sizeof(beamdata));
  if(obj==NULL) {
    fprintf(stderr,"primary_beam: no malloc for main struct\n");
    return NULL;
  }

   /* read-in Fourier-transformed beam file for this frequency channel and polarization */
	obj->real = create5DMatrixFloat(2,NUM_W_PLANES,NUM_POINT,SIZE_BEAM,SIZE_BEAM);
	obj->imag = create5DMatrixFloat(2,NUM_W_PLANES,NUM_POINT,SIZE_BEAM,SIZE_BEAM);
	/* No w-splitting settings
  obj->real = create3DMatrixFloat(NUM_POINT,SIZE_BEAM,SIZE_BEAM);
  obj->imag = create3DMatrixFloat(NUM_POINT,SIZE_BEAM,SIZE_BEAM);
	*/

  obj->u_tab = calloc(SIZE_BEAM,sizeof(float));
    obj->normal = create2DMatrixFloat(NUM_W_PLANES,NUM_POINT);
    obj->total_norm_sq = create2DMatrixFloat(NUM_W_PLANES,NUM_POINT);


    /* original RBW beam */
    
    //     if (pol == 0) sprintf(beamfilename,"%sbeam_ft_%03d_9_15_xx_sq.dat",getenv("BEAMDIR"),freq_index);
    //     if (pol == 1) sprintf(beamfilename,"%sbeam_ft_%03d_9_15_yy_sq.dat",getenv("BEAMDIR"),freq_index);
    
    /* New Curtin beam!! */
    
		    if (pol == 0) sprintf(beamfilename,"%sbeam_ft_%03d_9_15_xx_sq_new.dat",getenv("BEAMDIR"),freq_index);
		    if (pol == 1) sprintf(beamfilename,"%sbeam_ft_%03d_9_15_yy_sq_new.dat",getenv("BEAMDIR"),freq_index);
   
    
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
		for (k=0;k<NUM_W_PLANES;k++){	
		for (l=0;l<NUM_POINT;l++){

  		for (ii=0;ii<SIZE_BEAM;ii++){
    		for (j=0;j<SIZE_BEAM;j++){
      			if (fread(&obj->real[r][k][l][j][ii],sizeof(float),1,fp) != 1) printf("Error reading beam file.\n");
      			if (fread(&obj->imag[r][k][l][j][ii],sizeof(float),1,fp) != 1) printf("Error reading beam file.\n");

    		}
  		}

		}
		}
		}
		fclose(fp);

    for (k=0;k<NUM_W_PLANES;k++){
    for (l=0;l<NUM_POINT;l++){
    
        for (ii=0;ii<SIZE_BEAM;ii++){
            for (j=0;j<SIZE_BEAM;j++){
        
                obj->normal[k][l] += sqrt((obj->real[0][k][l][j][ii])*(obj->real[0][k][l][j][ii]) + (obj->imag[0][k][l][j][ii])*(obj->imag[0][k][l][j][ii]));
               // /(DELTA_U/INTRINSIC_DELTA_U)/(DELTA_U/INTRINSIC_DELTA_U);
                obj->total_norm_sq[k][l] += sqrt((obj->real[1][k][l][j][ii])*(obj->real[1][k][l][j][ii]) + (obj->imag[1][k][l][j][ii])*(obj->imag[1][k][l][j][ii]));
               //  /(DELTA_U/INTRINSIC_DELTA_U)/(DELTA_U/INTRINSIC_DELTA_U);
            
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
 ! NAME:		get_pb_values
 ! PURPOSE:		get primary beam values at the specified uv points in the grid
 ! ARGUMENTS:		pointer to delay array
 ! RETURNS:		integer (0=success)
******************************/

int get_pb_values(beamdata *bdata, int w_index, int point, int kernel_size, double *uvalues, double *vvalues, double uu, double vv, double *beam_real, double *beam_imag, float factor_scale) {

beam_diff_mwa_vector_2D(bdata->real[0][w_index][point],bdata->imag[0][w_index][point],bdata->u_tab,uvalues,uu,vvalues,vv,kernel_size,beam_real,beam_imag,factor_scale);

/* No w-splitting */
//beam_diff_mwa_vector_2D(bdata->real[point],bdata->imag[point],bdata->u_tab,uvalues,uu,vvalues,vv,kernel_size,beam_real,beam_imag);

return 0;


}


/******************************
 ! NAME:		get_pbsq_values
 ! PURPOSE:		get primary beam values at the specified uv points in the grid
 ! ARGUMENTS:		pointer to delay array
 ! RETURNS:		integer (0=success)
******************************/

int get_pbsq_values(beamdata *bdata, int w_index, int point, int kernel_size, double *uvalues, double *vvalues, double uu, double vv, double *beamsq_real, double *beamsq_imag, float factor_scale) {
    
beam_diff_mwa_vector_2D(bdata->real[1][w_index][point],bdata->imag[1][w_index][point],bdata->u_tab,uvalues,uu,vvalues,vv,kernel_size,beamsq_real,beamsq_imag,factor_scale);

/* No w-splitting */

//beam_diff_mwa_vector_2D(bdata->real[1][0][point],bdata->imag[1][0][point],bdata->u_tab,uvalues,uu,vvalues,vv,kernel_size,beamsq_real,beamsq_imag,factor_scale);

//beam_diff_mwa_vector_2D(bdata->real[point],bdata->imag[point],bdata->u_tab,uvalues,uu,vvalues,vv,kernel_size,beam_real,beam_imag);

return 0;


}

/******************************
 ! NAME:		free_pb
 ! PURPOSE:		frees the beam data structure
 ! ARGUMENTS:		beamdata structure type
 ! RETURNS:		void
******************************/


void free_pb(beamdata *b){

 int j,k,l,r;

if (b->u_tab !=NULL) free(b->u_tab);

for (r=0;r<2;r++){
for (k=0;k<NUM_W_PLANES;k++){
	for (l=0;l<NUM_POINT;l++){
		for (j=0;j<SIZE_BEAM;j++){
			free(b->real[r][k][l][j]);
			free(b->imag[r][k][l][j]);
	//		free(b->real[l][j]);
	//		free(b->imag[l][j]);

		}

	}
}
}
    
    for (k=0;k<NUM_W_PLANES;k++){
    free(b->normal[k]);
    free(b->total_norm_sq[k]);
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

  /* interpolate beam at given uprime,vprime values */
//    factor_scale=1.;
for (i=0;i<size;i++){
    utemp = u[i];///factor_scale;
    vtemp = v[i];///factor_scale;
	if ((sqrt((utemp-uprime)*(utemp-uprime)) >= u_tab[0])||(sqrt((vtemp-vprime)*(vtemp-vprime)))) {
	  xloc = ((utemp-uprime)*factor_scale-u_tab[0])/(u_tab[1]-u_tab[0]);
	  yloc = ((vtemp-vprime)*factor_scale-u_tab[0])/(u_tab[1]-u_tab[0]);
	   
	interp_bilin(beam_r,xloc,yloc,&out_real[i]);
	interp_bilin(beam_i,xloc,yloc,&out_imag[i]);

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
 ! NAME:		interp
 ! PURPOSE:		2D cublic spline interpolation for complex numbers
   sx,sy: size of grid x and y axes
   u,v: desired location to interpolate to in pixel units.

***************************************************/

#define SIZ_SPLINE 3

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

  /*  if (debug) {
    printf("interp_vis: p,q,m,n: %g,%g,%d,%d \n",p,q,m,n);
    for(j=0; j<SIZ_SPLINE; j++) {
      for(i=0; i<SIZ_SPLINE; i++) {
	val = getUVCell2D(i+m-1,j+n-1,sx,sy, uvgrid);
        printf("(%g,%g) ",creal(val),cimag(val));
      }
      printf("\n");
    }
  }
  */

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


    /* set y values for spline (REAL) */
    for (j=0; j<SIZ_SPLINE; j++){         /* loop over rows */
      for (i=0; i<SIZ_SPLINE; i++){       /* loop over cols */
	y[j][i]=uvgrid[m+i-1][n+j-1];
      }
    }

    /* initialise spline coeffs */
    splie2(x1,x2,y,SIZ_SPLINE,SIZ_SPLINE,y2);
    /* now do the interp (REAL) */
    splin2(x1,x2,y,y2,SIZ_SPLINE,SIZ_SPLINE,q,p,&re);


  *res_re = re;
    
    free(x1);
    free(x2);
    free(y);
    free(y2);
    
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
    
    *res_re = re;
    
    return 0;
}


/************************************************/

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



