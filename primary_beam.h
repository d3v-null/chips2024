/* structure definitions and utilities for primary beam calculation
Cathryn Trott. August 2013.
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "cspline.h"
//#include "suitesparse/cholmod.h"

#define VERSION 3.0
//#define SIZE_BEAM 151
#define SIZE_BEAM 95
#define NUM_FINE 16
#define SIZE_BEAM_FULL 285
//#define SIZE_BEAM 33
#define SIZE_SOURCE_NAME 16
#define SIZE_POL_TYPE 4
#define DELTA_U 0.5
#define DELTA_Us 0.5
#define DELTAUU 5.
#define SIZE_SIM 1200
//#define DELTA_U 1.0
#define INTRINSIC_DELTA_U 0.166
#define UMAX 300.
#define LOWER_FREQ_ULTRA 74.875e6
#define LOWER_FREQ 138.875e6
//#define LOWER_FREQ 115.e6
// FHD#define LOWER_FREQ 138.875e6
#define LOWER_FREQ_HIGH 167.035e6   //RTS
#define LOWER_FREQ_HIGH_EOR2 167.075e6   //RTS
//#define LOWER_FREQ_HIGH 167.075e6  //FHD
#define COARSE_CHAN_WIDTH 1.28e6
//#define DELTA_WSTACK 30.
#define WSTACK_DELT_TOT 2.
#define DELTA_WSTACK 2.
#define DELTA_W 0.
//#define DELTA_W 0.
#define HA_MAX 27.*M_PI/180.
#define LATITUDE -26.70331940
#define LONGITUDE 116.67081524
#define NPOL 2
//#define NUM_W_STACK 15
//#define NUM_W_PLANES 15
#define NUM_W_STACK 15
#define NUM_W_PLANES 15
#define WMAX ((NUM_W_PLANES/2)*DELTA_W + (NUM_W_STACK)*DELTA_WSTACK)
#define NUM_POINT 9
#define BEAM_SIZE_FLOAT 31.
//#define BEAM_SIZE_FLOAT 31.
#define HALF_PLANE_FLAG 0

#define CHAN_WIDTH 80.e3
#define PERIOD 8.
#define FIELD 0

extern float global_chanwidth;
extern float global_period;
extern int field;
extern float global_umax;
extern float global_lower_freq;

/* Micro w settings */
//#define NUM_W_PLANES 1
//#define NUM_W_STACK 20
//#define DELTA_W 4.
//#define DELTA_WSTACK 4.

/* structure that contains the original beam dataset and the tabulated u values */
typedef struct {
  float *****real;      /* array the size of NUM_W_PLANES x NUM_POINT x size_beam x size_beam FLOATS containing the real part of the beam */
  float *****imag;   /* array the size of NUM_W_PLANES x NUM_POINT x size_beam x size_beam FLOATS containing the imag part of the beam */
//  float ***real;      /* array the size of NUM_POINT x size_beam x size_beam FLOATS containing the real part of the beam */
//  float ***imag;   /* array the size of NUM_POINT x size_beam x size_beam FLOATS containing the imag part of the beam */
  float *u_tab;     /* array the size of size_beam FLOATS containing the tabulated u values for the beam */
    float **total_norm_sq;   /* float array containing the full normalised area under the squared beam */
    float **normal;   /* float array containing the full normalised area under the beam */
} beamdata;

/* structure that contains the original beam dataset and the tabulated u values */
typedef struct {
    float ******real;      /* array the size of NUM_W_PLANES x NUM_POINT x size_beam x size_beam FLOATS containing the real part of the beam */
    float ******imag;   /* array the size of NUM_W_PLANES x NUM_POINT x size_beam x size_beam FLOATS containing the imag part of the beam */
    //  float ***real;      /* array the size of NUM_POINT x size_beam x size_beam FLOATS containing the real part of the beam */
    //  float ***imag;   /* array the size of NUM_POINT x size_beam x size_beam FLOATS containing the imag part of the beam */
    float *u_tab;     /* array the size of size_beam FLOATS containing the tabulated u values for the beam */
    float ***total_norm_sq;   /* float array containing the full normalised area under the squared beam */
    float ***normal;   /* float array containing the full normalised area under the beam */
} beamdataf;


/* public function prototypes */
float ****create4DMatrixFloat(int sizex, int sizey, int sizez, int sizew);
float *****create5DMatrixFloat(int sizer, int sizex, int sizey, int sizez, int sizew);
float ******create6DMatrixFloat(int sizech, int sizer, int sizex, int sizey, int sizez, int sizew);
float ***create3DMatrixFloat(int sizex, int sizey, int sizez);
float **create2DMatrixFloat(int sizex, int sizey);
beamdataf *init_pbf(int freq_index, int pol);
beamdata *init_pb(int freq_index, int pol);
beamdata *init_pb_gauss(int freq_index, int pol);
void free_pb(beamdata *b);
void free_pbf(beamdataf *b);
int init_pb_delays(float delays);
int get_pb_values(beamdata *data, int w_index, int point, int kernel_size, double *uvalues, double *vvalues, double uu, double vv, double *beam_real, double *beam_imag, float factor_scale);
int get_pbsq_values(beamdata *data, int w_index, int point, int kernel_size, double *uvalues, double *vvalues, double uu, double vv, double *beamsq_real, double *beamsq_imag, float factor_scale);
int get_pbsq_values_fine(beamdataf *data, int w_index, int point, int kernel_size, double *uvalues, double *vvalues, double uu, double vv, double *beamsq_real, double *beamsq_imag, float factor_scale, int ch);
//int calc_bb_bv(uvdata *data, beamdata *beam_data, int size_fin, double WMAX, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, int flag_set);
//int calc_bb_bv_smallbeam(uvdata *data, beamdata *beam_data, int size_fin, double WMAX, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, int flag_set, unsigned long obsid);
int calc_bv_flags(uvdata *data, beamdata *beam_data, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, int flag_set, unsigned long obsid, char *ext);
int calc_bv_flags_uv_fine(uvdata *data, uvdata *data2, beamdataf *beam_data, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext, int *flag14m);
int calc_bv_flags_uv_sims(uvdata *data, uvdata *data2, float **data_sim, float **data_simi, beamdataf *beam_data, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext, int *flag14m, int coarsenum);
int calc_bv_flags_diff(uvdata *data, uvdata *data2, beamdata *beam_data, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext, int *flag14m);
int calc_bv_gauss_analytic(uvdata *data, uvdata *data2, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext, int *flag14m);
int calc_bv_bn_analytic(uvdata *data, uvdata *data2, int size_fin, int ch, double frequency, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext, int *flag14m, int band, float umx);
int calc_bv_bn_stack(uvdata *data, uvdata *data2, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int pol, int u_size, double systemp, unsigned long obsid, char *ext, int *flag14m, int band);
int calc_bv_bn_simmo(uvdata *data, uvdata *data2, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int pol, int u_size, double systemp, unsigned long obsid, char *ext, int *flag14m, int band);
int calc_bv_flags_gauss(uvdata *data, uvdata *data2, beamdata *beam_data, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext, int *flag14m);
int calc_bv_flags_diff_fine(uvdata *data, uvdata *data2, beamdataf *beam_data, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext, int *flag14m);
int calc_bv_flags_ultra(uvdata *data, uvdata *data2, beamdata *beam_data, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext);
int calc_bv_flags_check(uvdata *data, uvdata *data2, beamdata *beam_data, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext);
int calc_bv_sims(beamdata *beam_data, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, char *ext, int fg);
int calc_bv_flags_sims(uvdata *data, uvdata *data2, beamdata *beam_data, int size_fin, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, unsigned long obsid, char *ext, float S1, float S2);
int calc_grid(double **uv, int size_fin, int ch, double frequency, int point, int freq_index, int wstack, int pol, int u_size, int nbaselines, float bp, int ext);
//float complex *check_cal_func(uvdata *data, int pol, int updown);
//int make_flags(uvdata *data, beamdata *beam_data, int size_fin, double WMAX, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, int flag_set, unsigned long obsid);
//int calc_bb_bv_grid(uvdata *data, int size_fin, double WMAX, int ch, double l_rot, double m_rot, double frequency, int point, int freq_index, int wstack, int pol, int u_size, double systemp, int flag_set);

