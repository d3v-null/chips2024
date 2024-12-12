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
#include "uvfits.h"
#include "primary_beam.h"
#include "fitsio.h"
#include "complex.h"

// gcc -Wall -g grid_vis.c cspline.c uvfits.c -o gridvis -L/usr/local/lib -lm -lcfitsio -lcholmod

/* Function prototypes */
double gauss_noise2(double value1);
void free_mem3D(double ***matrix, int xsize, int ysize, int zsize);
double vectorVectorMultiply(int vsize, double vec1[], double vec2[]);
float complex **create2DMatrixFloatComplex(int sizex, int sizey);
void free_memFloat(float **matrix, int xsize, int ysize);
int matrixVectorMultiply(int vsize, double vec[], double **mat, double outVector[]);
double **create2DMatrix(int sizex, int sizey);
double ***create3DMatrix(int sizex, int sizey, int sizez);
int ***create3DMatrixInt(int sizex, int sizey, int sizez);
long ***create3DMatrixLong(int sizex, int sizey, int sizez);
double ****create4DMatrix(int sizex, int sizey, int sizez, int sizew);
double complex **create2DMatrixComplex(int sizex, int sizey);
double complex ***create3DMatrixComplex(int sizex, int sizey, int sizez);
void memzero2DFloat(float **matrix, int xsize, int ysize);
long ****create4DMatrixLong(int sizex, int sizey, int sizez, int sizew);
// float **create2DMatrixFloat(int sizex, int sizey);
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
void memzero2D(double **matrix, int xsize, int ysize);
void memzero2Dint(int **matrix, int xsize, int ysize);
void memzero2Dlong(long int **matrix, int xsize, int ysize);
void memzero3D(double ***matrix, int xsize, int ysize, int zsize);
void memzero3Dint(int ***matrix, int xsize, int ysize, int zsize);
void free_mem(double **matrix, int xsize, int ysize);

/* global variables */
float global_period = PERIOD;
float global_chanwidth = CHAN_WIDTH;
int field = FIELD;
int debug = 0;
// char *infilestub=NULL;
char *infilelist = NULL;
FILE *fpd;

void usage()
{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "test_readuvfits <options> -i inputfilename\n");
  fprintf(stderr, "\t-d debug_level (0=no debugging)\n");
  exit(0);
}

void parse_cmdline(const int argc, char *const argv[])
{
  int result = 0;
  const char optstring[] = "d:i:p:c:f:bandlowchan";

  while ((result = getopt(argc, argv, optstring)) != -1)
  {
    switch (result)
    {
    case 'i':
      infilelist = optarg;
      break;
    case 'p':
      global_period = atof(optarg);
      break;
    case 'c':
      global_chanwidth = atof(optarg);
      break;
    case 'f':
      field = atoi(optarg);
      break;
    case 'd':
      debug = atoi(optarg);
      break;
    default:
      fprintf(stderr, "unknown option: %c\n", result);
      usage(argv[0]);
    }
  }
}

/* Main */
int main(int argc, char *argv[])
{

  int freq_index, ch, size_fin, flag_tsys = 0, i_coarse, i, j, k, flag, Nchunk = 8;
  FILE *fpobs = NULL, *ftsysobs = NULL, *flog = NULL, *fmask = NULL, *flut = NULL, *fptrv = NULL, *fptrflags = NULL, *fptrdiff = NULL, *fptruvw = NULL, *fptr = NULL;
  int pol, polindex, polarray[2] = {0, 1}, Nchan = 384, Neta = 192, nchansing = 16;
  double frequency, uu, vv, ww, norm, fac, v11, v22, rsq;
  double tsys[2];
  double ra_app, dec_app, ra_point, dec_point, ra_offset, dec_offset, ra_field, dec_field, ra_zenith, dec_zenith;
  double u_size, systemp, weight, v1, v2, v3, v4, numb2 = 0., numb = 0., numb3 = 0., denom1 = 1., denom2 = 1., denom3 = 1.;
  char syslogfile[1024], infilename[1024], maskfile[1024], filename_real[1024], filename_real2[1024], filename_flags[1024], debugfile[1024], filename_real3[1024], filename_uvw[1024], filename_noise[1024], infilestub[1024], *ext = NULL;
  char *obsid;
  char outfilename[1024], temp[1024], outfilename2[1024], outfilename3[1024];
  float vis_rtot, vis_itot, vis_rdiff, vis_idiff, **flag_array = NULL, *mask = NULL, az, el, dec_rad, mjd, lmst, m_rot, l_rot, az_point, alt_point, clight = 2.997e8;
  float beamsigma = 0.1, D = 4.4, vis_corr = 1., lowfreq = LOWER_FREQ, factor_kpa = 0., DM = 0., BW = 0., power_window = 0., power_wedge = 0., power_all = 0., **numdiff = NULL, **diff = NULL;
  int wstack = 0, point = 4, nthreads = 1, band = 0, nbase = 8128;
  float **uvw = NULL, *nut_func = NULL, ha_rad, ha_point, distance = 0., w1 = 0., w2 = 0., norm_factor = 1., ***visall = NULL;
  float *mean = NULL, *mn = NULL, *stddev = NULL, *numstd = NULL, **final_stdnorm = NULL, countoutlier = 0.;
  float complex **fourier = NULL, **vis_tot = NULL, **vis_actual = NULL, **vis_diff = NULL;
  float a_nut[4] = {0.3635819, 0.4891775, 0.1365995, 0.0106411};
  // float a_nut[4]={1.,0.,0.,0.};

  //  fitsfile *fptr,*fptri;
  int res = 0, chunk = 0;
  uvdata *data, *data2;
  uvReadContext *iter, *iter2;

  fpd = stderr;

  // Information

  if (argc < 2)
  {
    fprintf(stderr, "Usage: %s <options> obsid band Nchan subdirectory/uvfits_name_stub \n", argv[0]);
    fprintf(stderr, "\t -p period\n");
    fprintf(stderr, "\t -c chanwidth\n");
    fprintf(stderr, "\t -field fieldnum\n");
    exit(1);
  }

  printf("Program to calculate the beam weights for an MWA observation.\n\n");

  // infilestub = argv[1];
  //    infilelist = argv[1];
  obsid = (argv[1]);
  band = atoi(argv[2]);
  Nchan = atoi(argv[3]);
  ext = argv[4];

  parse_cmdline(argc, argv);

  {

    char tmp1[FILENAME_MAX], *p;

    /* Open log file */

    sprintf(syslogfile, "%ssyslog_checksum.txt", getenv("OUTPUTDIR"));
    //	printf("syslogfile: %s\n",syslogfile);
    if ((flog = fopen(syslogfile, "w")) == NULL)
    {
      fprintf(stderr, "Cannot open output log file\n");
      return 1;
    }
  }

  time_t tim = time(NULL);
  struct tm *tmptr = gmtime(&tim);

  fprintf(flog, "Processed with version %f of DS on %s\n", VERSION, asctime(localtime(&tim)));

  az = 0.;
  el = 90.0;

  /* define uv grid size and spacing */

  u_size = floor(UMAX / DELTA_U);

  if (HALF_PLANE_FLAG)
  {
    size_fin = (int)2. * u_size * u_size;
  }
  else
  {
    size_fin = (int)4. * u_size * u_size;
  }

  // printf("Num. u bins: %lg, size of uv plane: %d\n",u_size,size_fin);

  /* Build Fourier Kernel for LOS transform */

  Neta = Nchan / 2;
  float Nchanfloat = (float)Nchan;

  fourier = create2DMatrixFloatComplex(Neta, Nchan);
  for (j = 0; j < Nchan; j++)
  {
    //        for (i=0;i<Neta;i++) fourier[i][j] = cos(2.*M_PI*((float) i)*((float) j)/((float) Nchan)) - I*sin(2.*M_PI*((float) i)*((float) j)/((float) Nchan));
    for (i = 0; i < Neta; i++)
    {
      fourier[i][j] = cexp(-I * 2. * M_PI * ((float)i) * ((float)j) / (Nchanfloat)) / (float)Nchan;
      // printf("i %d j %d fourier %f %f\n",i,j,creal(fourier[i][j]),cimag(fourier[i][j]));
    }
  }

  nut_func = calloc(Nchan, sizeof(float));
  for (j = 0; j < Nchan; j++)
  {
    for (i = 0; i < 4; i++)
      nut_func[j] += a_nut[i] * cos(2. * M_PI * (float)i * (j - Nchan / 2.) / (float)Nchan);
  }

  //    printf("nut func %f %f\n",nut_func[0],nut_func[Nchan/2]);

  /* Define dense vector for bdag_v accumulation and flags for uv-sampling */

  nbase = 50000;

  /******************************************

  Begin loop over coarse channels

  *******************************************/

  sprintf(infilename, "%s/%s.uvfits", getenv("DATADIR"), ext);

  printf("Filename: %s\n", infilename);

  chunk = 0;

  /* read-in uvfits data */

  uvfitsSetDebugLevel(0);
  res = readUVFITSInitIterator(infilename, &data, &iter);
  res = readUVFITSInitIterator(infilename, &data2, &iter2);

  printf("%d\n", nbase);

  vis_tot = create2DMatrixFloatComplex(nbase, Neta);
  vis_actual = create2DMatrixFloatComplex(Nchan, 40);
  vis_diff = create2DMatrixFloatComplex(nbase, Neta);
  //   mean = create2DMatrixFloat(Nchan,Nchunk);
  //   stddev = create2DMatrixFloat(Nchan,Nchunk);
  //    mn = create2DMatrixFloat(Nchan,Nchunk);
  //   numstd = create2DMatrixFloat(Nchan,Nchunk);
  visall = create3DMatrixFloat(Nchan, Nchunk, nbase);
  diff = create2DMatrixFloat(nbase, Nchan);
  numdiff = create2DMatrixFloat(Nchan, Nchunk);
  final_stdnorm = create2DMatrixFloat(Nchan, Nchunk);

  uvw = create2DMatrixFloat(nbase, 3);

  /* Key parameters for cosmological transformation - magic numbers that are secrets of the cosmos */

  BW = 30.72e6;
  factor_kpa = 1.78496e-06;
  DM = 6200. / 2. / M_PI;

  for (i = 0; i < nbase; i++)
  {
    for (k = 0; k < Neta; k++)
    {
      vis_tot[i][k] = 0. + I * 0.;
      vis_diff[i][k] = 0. + I * 0.;
    }
  }

  if (res != 0)
  {
    fprintf(stderr, "readUVFITSInitIterator failed on timestep2 with error %d\n", res);
    return res;
  }

  fflush(stdout);

  /* Dummy read of data2 to increment it to second timestep */
  if ((res = readUVFITSnextIter(data2, iter2)) != 0)
  {
    //   printf("Second set iteration not found\n");
    return res;
  }

  /* code runs entirely within WHILE loop, which loops over sections of the input data file (one timestep per iteration) */

  uvfitsSetDebugLevel(0);
  printf("freq. channels %d, num vis %d\n", data->n_freq, data->n_vis);

  while ((res = readUVFITSnextIter(data, iter)) == 0)
  {
    printf("Chunk %d. Time: %f. baselines: %d\n", chunk++, data->date[0], data->n_baselines[0]);

    if (band == 0)
      lowfreq = LOWER_FREQ;
    if (band == 1)
      lowfreq = LOWER_FREQ_HIGH;

    /************************************************/
    /* Loop over channels to compute each frequency separately for each OpenMP thread */

    for (i = 0; i < nbase; i++)
    {
      for (k = 0; k < Neta; k++)
      {
        vis_tot[i][k] = 0. + I * 0.;
        vis_diff[i][k] = 0. + I * 0.;
      }
    }

    for (ch = 0; ch < data->n_freq; ch++)
    {

      // get frequency information //

      freq_index = (((data->cent_freq) - data->n_freq / 2 * data->freq_delta) + ch * data->freq_delta - LOWER_FREQ) / COARSE_CHAN_WIDTH;
      frequency = ((data->cent_freq) - data->n_freq / 2 * data->freq_delta) + ch * data->freq_delta;

      //	printf("Frequency: %g Channel: %d %f %f\n",frequency,ch,CHAN_WIDTH,data->freq_delta);

      pol = 0;

      for (i = 0; i < data->n_baselines[0]; i++)
      {
        //         printf("baseline: %d\n",i);

        //       if (data->weightdata[0][(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)] > 64.) printf("weight large!\n");

        // contrib_vis++;

        uu = data->u[0][i] * (frequency);
        vv = data->v[0][i] * (frequency);
        ww = data->w[0][i] * (frequency);
        flag = 0;

        distance = sqrt(uu * uu + vv * vv);
        // printf("Freq: %g, Distance: %g\n",frequency,distance);

        // add and difference two visibility sets in real and imaginary

        vis_itot = 0.;
        vis_rtot = 0.;
        vis_idiff = 0.;
        vis_rdiff = 0.;

        // Compute the normalising factor to handle different weights

        w1 = data->weightdata[0][(i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol)];
        w2 = data2->weightdata[0][(i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol)];

        norm_factor = sqrt(1. / w1 + 1. / w2);

        //      if ((ch == 380)&&(i == 99)&&(pol == 1)) printf("Is it the same? %f %f\n",(data->visdata[0][2*(i*(data->n_pol*data->n_freq) + ch*data->n_pol + pol)]),(data2->visdata[0][2*(i*(data2->n_pol*data2->n_freq) + ch*data2->n_pol + pol)]));

        if ((data->weightdata[0][(i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol)] > 0.) && (data2->weightdata[0][(i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol)] > 0.))
          vis_rtot = (data->visdata[0][2 * (i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol)] + data2->visdata[0][2 * (i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol)]);

        if ((data->weightdata[0][(i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol + 1)] > 0.) && (data2->weightdata[0][(i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol + 1)] > 0.))
          vis_itot = (data->visdata[0][2 * (i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol) + 1] + data2->visdata[0][2 * (i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol) + 1]);

        if ((data->weightdata[0][(i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol)] > 0.) && (data2->weightdata[0][(i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol)] > 0.))
          vis_rdiff = (data->visdata[0][2 * (i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol)] - data2->visdata[0][2 * (i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol)]);

        if ((data->weightdata[0][(i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol + 1)] > 0.) && (data2->weightdata[0][(i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol + 1)] > 0.))
          vis_idiff = (data->visdata[0][2 * (i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol) + 1] - data2->visdata[0][2 * (i * (data->n_pol * data->n_freq) + ch * data->n_pol + pol) + 1]);

        if (debug)
          fflush(flog);

        /* compute u and v locations for each cell */

        /* Loop over output LOS Fourier bins and contribute to each */
        for (k = 0; k < Neta; k++)
        { // loop over eta modes

          vis_tot[i][k] += fourier[k][ch] * (vis_rtot + I * vis_itot) * nut_func[ch]; //*vis_corr;
          vis_diff[i][k] += fourier[k][ch] * (vis_rdiff + I * vis_idiff) * nut_func[ch];

          // if ((i == 0)&&(k == 0)) printf("Outputs: %f %f %f %f\n",fourier[k][ch],vis_rtot,vis_itot,creal(vis_tot[i][k]));
        }
      }

      //     printf("Channel %d Chunk %d\n",ch,chunk);

      // printf("Number of cells: %d %f %f %f\n",i,numb,numb2,distance);

      fflush(flog);
      //      printf("Herech\n");
    } //  ******** end loop over channels ********
      //      printf("Here1\n");
    chunk++;

    for (i = 0; i < data->n_baselines[0]; i++)
    {

      uu = data->u[0][i] * (frequency);
      vv = data->v[0][i] * (frequency);
      ww = data->w[0][i] * (frequency);

      distance = sqrt(uu * uu + vv * vv);

      for (k = 0; k < Neta; k++)
      { // loop over eta modes

        if (((float)k > factor_kpa * BW / DM * distance) && (k < 22) && (distance > 0) && (distance < 100))
        {
          //		printf("1 %d %f %f\n",k,distance,factor_kpa*BW/DM*distance);
          power_window += vis_tot[i][k] * conj(vis_tot[i][k]) - vis_diff[i][k] * conj(vis_diff[i][k]);
          numb += 1.;
        }
        if (((float)k < factor_kpa * BW / DM * distance) && (distance > 0) && (distance < 100))
        {
          //		printf("2 %d %f %f\n",k,distance,factor_kpa*BW/DM*distance);
          power_wedge += vis_tot[i][k] * conj(vis_tot[i][k]) - vis_diff[i][k] * conj(vis_diff[i][k]);
          numb2 += 1.;
        }

        power_all += vis_tot[i][k] * conj(vis_tot[i][k]) - vis_diff[i][k] * conj(vis_diff[i][k]);
        numb3 += 1.;
      }
    }

    if ((res = readUVFITSnextIter(data2, iter2)) != 0)
    {
      printf("Second set iteration not found\n");
      res = 0;
      ffcmrk;
      break;
    }
    if ((res = readUVFITSnextIter(data, iter)) != 0)
    {
      printf("Second set iteration not found\n");
      res = 0;
      ffcmrk;
      break;
    }

    //     printf("Here2\n");

  } //*************** END WHILE LOOP OVER UVFITS TIME STEPS *******/

  if (res != 1)
  {
    fprintf(stderr, "readUVFISTnextIter returned %d\n", res);
  }

  printf("Iter closing\n");

  readUVFITSCloseIter(iter);

  // printf("Iter2 closing\n");
  readUVFITSCloseIter(iter2);
  // printf("Iter/Iter2 closed\n");

  /* Output FT-ed visibilities */

  /************** OUTPUT FILES ********************/
  /* define filenames */

  denom1 = numb;
  denom2 = numb2;
  denom3 = numb3;

  // printf("Number of cells: %f %f\n",numb,numb2);

  sprintf(filename_real2, "%soutput_metrics_%s.dat", getenv("OUTPUTDIR"), obsid);
  printf("%s\n", filename_real2);

  // Check to see if the file already exists

  if ((fptrv = fopen(filename_real2, "a+")) != NULL)
  {

    fprintf(fptrv, "%s %g %g %g %g %g %g\n", obsid, power_wedge / denom2, denom2, power_window / denom1, denom1, power_all / denom3, denom3);

    fclose(fptrv);
  }

  power_wedge = 0.;
  power_window = 0.;
  power_all = 0.;
  numb = 0.;
  numb2 = 0.;
  numb3 = 0.;

  free(data);
  free(data2);

  fclose(flog);

  return 0;
}

// Functions

double gauss_noise2(/* value1, value2) */
                    double value1)
{

  double fac, rsq, v1, v2;

  while (1)
  {
    // Random location in unit square
    v1 = 2.0 * drand48() - 1.0;
    v2 = 2.0 * drand48() - 1.0;
    rsq = v1 * v1 + v2 * v2;
    // Take only those inside unit circle
    if (rsq < 1.0)
      break;
  }
  // Following the recipe ...
  fac = sqrt(-2.0 * log(rsq) / rsq);
  //    printf("gauss_noise v1 %lg rsq %lg fac %lg\n",v1,rsq,fac);
  return v1 * fac;
}

double trace_mat(double **matrix, int size)
{
  int i;
  double out = 0.;

  for (i = 0; i < size; i++)
    out += matrix[i][i];

  //  printf("out %g\n",out);
  return out;
}

/*
int round(float num){

int n = (int)(num < 0 ? (num - 0.5) : (num + 0.5));

 return n;
}
*/

void where(long *vec_in, long int size, int limit, long int *vec_out, int *counter)
{
  int i, c = 0;

  for (i = 0; i < size; i++)
  {
    if (vec_in[i] == limit)
    {
      vec_out[c] = i;
      c++;
    }
  }

  counter = &c;
  //     printf("c %i\n",c);
}

int matrixVectorMultiply(int vsize, double vec[], double **mat, double outVector[])
{
  int i, j;

  for (i = 0; i < vsize; i++)
  {
    outVector[i] = 0.;
    for (j = 0; j < vsize; j++)
    {
      outVector[i] += vec[j] * mat[i][j];
    }
  }
  return 0;
}

int matrixMatrixMultiply(int x1size, int y1size, int x2size, int y2size, double **mat1, double **mat2, double **outMat)
{
  int i, j, k;

  if (x1size != y2size)
    exit(1);

  for (i = 0; i < x2size; i++)
  {
    for (j = 0; j < y1size; j++)
    {
      for (k = 0; k < x1size; k++)
      {

        outMat[i][j] += mat1[k][j] * mat2[i][k];
      }
    }
  }
  return 0;
}

int CmatrixMatrixMultiply(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag)
{
  int i, j, k;

  if (x1size != y2size)
    exit(1);

  for (i = 0; i < x2size; i++)
  {
    for (j = 0; j < y1size; j++)
    {
      for (k = 0; k < x1size; k++)
      {

        outMatreal[i][j] += mat1real[k][j] * mat2real[i][k] - mat1imag[k][j] * mat2imag[i][k];
        outMatimag[i][j] += mat1real[k][j] * mat2imag[i][k] + mat1imag[k][j] * mat2real[i][k];
      }
    }
  }
  return 0;
}

int CmatrixMatrixMultiplyConj(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag)
{
  int i, j, k;

  /* Conjugates the imaginary part of the second matrix */

  if (x1size != y2size)
    exit(1);

  for (i = 0; i < x2size; i++)
  {
    for (j = 0; j < y1size; j++)
    {
      for (k = 0; k < x1size; k++)
      {

        outMatreal[i][j] += mat1real[k][j] * mat2real[i][k] + mat1imag[k][j] * mat2imag[i][k];
        outMatimag[i][j] += -mat1real[k][j] * mat2imag[i][k] + mat1imag[k][j] * mat2real[i][k];
      }
    }
  }
  return 0;
}

int CmatrixMatrixMultiplyConjTrans(int x1size, int y1size, int x2size, int y2size, double **mat1real, double **mat1imag, double **mat2real, double **mat2imag, double **outMatreal, double **outMatimag)
{
  int i, j, k;

  /* Conjugates the imaginary part of the second matrix, and transposes */

  /* Must be symmetric!! */

  if (x1size != y2size)
    exit(1);

  for (i = 0; i < x2size; i++)
  {
    for (j = 0; j < y1size; j++)
    {
      for (k = 0; k < x1size; k++)
      {

        outMatreal[i][j] += mat1real[k][j] * mat2real[k][i] + mat1imag[k][j] * mat2imag[k][i];
        outMatimag[i][j] += -mat1real[k][j] * mat2imag[k][i] + mat1imag[k][j] * mat2real[k][i];
      }
    }
  }
  return 0;
}

int matrixTranspose(int xsize, int ysize, double **mat_in, double **mat_out)
{
  int i, j;

  for (i = 0; i < xsize; i++)
  {
    for (j = 0; j < ysize; j++)
    {

      mat_out[j][i] = mat_in[i][j];
    }
  }
  return 0;
}

double vectorVectorMultiply(int vsize, double vec1[], double vec2[])
{
  int i;
  double temp;

  temp = 0.;
  for (i = 0; i < vsize; i++)
  {
    temp += vec1[i] * vec2[i];
  }
  return temp;
}

double **create2DMatrix(int sizex, int sizey)
{
  double **output;
  int k;

  output = calloc(sizex, sizeof(double *));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(double));
  }
  return output;
}

double ***create3DMatrix(int sizex, int sizey, int sizez)
{
  double ***output;
  int k, l;

  output = calloc(sizex, sizeof(double **));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(double *));
  }
  for (k = 0; k < sizex; k++)
  {
    for (l = 0; l < sizey; l++)
    {
      output[k][l] = calloc(sizez, sizeof(double));
    }
  }
  return output;
}

int ***create3DMatrixInt(int sizex, int sizey, int sizez)
{
  int ***output;
  int k, l;

  output = calloc(sizex, sizeof(int **));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(int *));
  }
  for (k = 0; k < sizex; k++)
  {
    for (l = 0; l < sizey; l++)
    {
      output[k][l] = calloc(sizez, sizeof(int));
    }
  }
  return output;
}

long ***create3DMatrixLong(int sizex, int sizey, int sizez)
{
  long ***output;
  int k, l;

  output = calloc(sizex, sizeof(long **));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(long *));
  }
  for (k = 0; k < sizex; k++)
  {
    for (l = 0; l < sizey; l++)
    {
      output[k][l] = calloc(sizez, sizeof(long));
    }
  }
  return output;
}

double ****create4DMatrix(int sizex, int sizey, int sizez, int sizew)
{
  double ****output;
  int k, l, m;

  output = calloc(sizex, sizeof(double ***));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(double **));
  }
  for (k = 0; k < sizex; k++)
  {
    for (l = 0; l < sizey; l++)
    {
      output[k][l] = calloc(sizez, sizeof(double *));
    }
  }
  for (k = 0; k < sizex; k++)
  {
    for (l = 0; l < sizey; l++)
    {
      for (m = 0; m < sizez; m++)
      {
        output[k][l][m] = calloc(sizew, sizeof(double));
      }
    }
  }
  return output;
}

long ****create4DMatrixLong(int sizex, int sizey, int sizez, int sizew)
{
  long ****output;
  int k, l, m;

  output = calloc(sizex, sizeof(long ***));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(long **));
  }
  for (k = 0; k < sizex; k++)
  {
    for (l = 0; l < sizey; l++)
    {
      output[k][l] = calloc(sizez, sizeof(long *));
    }
  }
  for (k = 0; k < sizex; k++)
  {
    for (l = 0; l < sizey; l++)
    {
      for (m = 0; m < sizez; m++)
      {
        output[k][l][m] = calloc(sizew, sizeof(long));
      }
    }
  }
  return output;
}

float **create2DMatrixFloat(int sizex, int sizey)
{
  float **output;
  int k;

  output = calloc(sizex, sizeof(float *));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(float));
  }
  return output;
}

/**************************************
 ! NAME:		create4DMatrix
 ! PURPOSE:		create a 4D double array
 ! RETURNS:		pointer to a 4D Float array
**************************************/

float ****create4DMatrixFloat(int sizex, int sizey, int sizez, int sizew)
{
  float ****output;
  int k, l, m;

  output = calloc(sizex, sizeof(float ***));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(float **));
  }
  for (k = 0; k < sizex; k++)
  {
    for (l = 0; l < sizey; l++)
    {
      output[k][l] = calloc(sizez, sizeof(float *));
    }
  }
  for (k = 0; k < sizex; k++)
  {
    for (l = 0; l < sizey; l++)
    {
      for (m = 0; m < sizez; m++)
      {
        output[k][l][m] = calloc(sizew, sizeof(float));
      }
    }
  }
  return output;
}

/**************************************
 ! NAME:		create3DMatrix
 ! PURPOSE:		create a 3D double array
 ! RETURNS:		pointer to a 3D Float array
**************************************/

float ***create3DMatrixFloat(int sizex, int sizey, int sizez)
{
  float ***output;
  int k, l;

  output = calloc(sizex, sizeof(float **));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(float *));
  }
  for (k = 0; k < sizex; k++)
  {
    for (l = 0; l < sizey; l++)
    {
      output[k][l] = calloc(sizez, sizeof(float));
    }
  }

  return output;
}

/************************************************/

double complex **create2DMatrixComplex(int sizex, int sizey)
{
  double complex **output;
  int k;

  output = calloc(sizex, sizeof(double complex *));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(double complex));
  }
  return output;
}

float complex **create2DMatrixFloatComplex(int sizex, int sizey)
{
  float complex **output;
  int k;

  output = calloc(sizex, sizeof(float complex *));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(float complex));
  }
  return output;
}

int **create2DMatrixInt(int sizex, int sizey)
{
  int **output;
  int k;

  output = calloc(sizex, sizeof(int *));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(int));
  }
  return output;
}

long int **create2DMatrixLong(int sizex, int sizey)
{
  long int **output;
  int k;

  output = calloc(sizex, sizeof(long int *));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(long int));
  }
  return output;
}

unsigned long **create2DMatrixShort(int sizex, int sizey)
{
  unsigned long **output;
  int k;

  output = calloc(sizex, sizeof(unsigned long *));
  for (k = 0; k < sizex; k++)
  {
    output[k] = calloc(sizey, sizeof(unsigned long));
  }
  return output;
}

double *create1DVector(int size)
{
  double *output;

  output = calloc(size, sizeof(double));
  return output;
}

void memzero2D(double **matrix, int xsize, int ysize)
{
  int i;

  for (i = 0; i < xsize; i++)
    memset(matrix[i], 0., ysize * sizeof(double));
}

void memzero2DFloat(float **matrix, int xsize, int ysize)
{
  int i;

  for (i = 0; i < xsize; i++)
    memset(matrix[i], 0., ysize * sizeof(float));
}

void memzero2Dint(int **matrix, int xsize, int ysize)
{
  int i;

  for (i = 0; i < xsize; i++)
    memset(matrix[i], 0., ysize * sizeof(int));
}

void memzero2Dlong(long int **matrix, int xsize, int ysize)
{
  int i;

  for (i = 0; i < xsize; i++)
    memset(matrix[i], 0., ysize * sizeof(long int));
}

void memzero3D(double ***matrix, int xsize, int ysize, int zsize)
{
  int i, j;

  for (i = 0; i < xsize; i++)
  {
    for (j = 0; j < ysize; j++)
      memset(matrix[i][j], 0., zsize * sizeof(double));
  }
}

void memzero3Dint(int ***matrix, int xsize, int ysize, int zsize)
{
  int i, j;

  for (i = 0; i < xsize; i++)
  {
    for (j = 0; j < ysize; j++)
      memset(matrix[i][j], 0., zsize * sizeof(int));
  }
}

void free_mem(double **matrix, int xsize, int ysize)
{
  int i;

  for (i = 0; i < xsize; i++)
    free(matrix[i]);
}

void free_memFloat(float **matrix, int xsize, int ysize)
{
  int i;

  for (i = 0; i < xsize; i++)
    free(matrix[i]);
}

void free_mem3D(double ***matrix, int xsize, int ysize, int zsize)
{
  int i, j;

  for (i = 0; i < xsize; i++)
  {
    for (j = 0; j < ysize; j++)
      free(matrix[i][j]);
  }
}
