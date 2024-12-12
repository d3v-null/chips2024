#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
// #include <complex.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>
#include <assert.h>
#include <libgen.h> // for basename()
#include "cspline.h"
#include <getopt.h>
#include "uvfits.h"
#include "primary_beam.h"
#include "fitsio.h"
#include <omp.h>

// gcc -Wall -g combine_data.c cspline.c uvfits.c -o combine_data -L/usr/local/lib -lm -lcfitsio

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
void memzero2DFloat(float **matrix, int xsize, int ysize);
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
void memzero2D(double **matrix, int xsize, int ysize);
void memzero2Dint(int **matrix, int xsize, int ysize);
void memzero2Dlong(long int **matrix, int xsize, int ysize);
void memzero3D(double ***matrix, int xsize, int ysize, int zsize);
void memzero3Dint(int ***matrix, int xsize, int ysize, int zsize);
void free_mem(double **matrix, int xsize, int ysize);

/* global variables */
int debug = 0;

/* global variables */
float global_period = PERIOD;
float global_chanwidth = CHAN_WIDTH;

char *infilename = NULL;
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
	const char optstring[] = "d:i:";

	while ((result = getopt(argc, argv, optstring)) != -1)
	{
		switch (result)
		{
		case 'i':
			infilename = optarg;
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

	int i, j, ch, size_fin, Nchan, size_half, band, addsub = 0, filenum = 0;
	FILE *fptrv = NULL, *fpbb = NULL, *flog = NULL, *fin = NULL;
	float frequency = LOWER_FREQ, *freq = NULL, *weights = NULL, *w1 = NULL;
	double chanwidth, freq_start;
	double u_size;
	char *date = NULL, fileext[100], syslogfile[1024], *pol = NULL, bvfilename[1024], bbfilename[1024], normfilename[1024], filename_real[1024], filename_real2[1024], filename_real6[1024], inputfile[1024], *ext = NULL;
	char weightsfilename[1024], vistfilename[1024], visdfilename[1024], noisetfilename[1024], noisedfilename[1024];
	unsigned long counter = 0;
	float complex *visdiff = NULL, *vistot = NULL, *noisediff = NULL, *noisetot = NULL, *visd = NULL, *vist = NULL, *noised = NULL, *noiset = NULL;

	fpd = stderr;

	// Information

	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s input_list Nchan ext add-0/subtract-1\n", argv[0]);
		exit(1);
	}

	printf("Program to combine data over multiple gridded sets.\n\n");

	date = argv[1];
	Nchan = atoi(argv[2]);
	ext = argv[3];
	addsub = atoi(argv[4]);

	/* Open log file */
	sprintf(syslogfile, "%ssyslog_combinedata_%s.txt", getenv("OUTPUTDIR"), ext);
	printf("syslogfile: %s\n", syslogfile);
	if ((flog = fopen(syslogfile, "w")) == NULL)
	{
		fprintf(stderr, "Cannot open output log file\n");
		return 1;
	}

	time_t tim = time(NULL);
	struct tm *tmptr = gmtime(&tim);

	fprintf(flog, "Processed with version %f of CHIPS on %s\n", VERSION, asctime(localtime(&tim)));

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

	size_half = (int)2. * u_size * (u_size + 1);
	size_fin = size_half;
	fprintf(flog, "size_half %d u_size %g Nchan %d \n", size_half, u_size, Nchan);

	weights = calloc(size_half * Nchan, sizeof(float));
	vistot = calloc(size_half * Nchan, sizeof(float complex));
	visdiff = calloc(size_half * Nchan, sizeof(float complex));
	noisetot = calloc(size_half * Nchan, sizeof(float complex));
	noisediff = calloc(size_half * Nchan, sizeof(float complex));

	filenum = 0;

	/****************************************************/

	// Open filename of input data and read
	sprintf(inputfile, "%s", date);
	printf("input file: %s\n", inputfile);
	if ((fin = fopen(inputfile, "r")) == NULL)
	{
		fprintf(stderr, "Cannot open input file\n");
		return 1;
	}

	while (fscanf(fin, "%s", fileext) != EOF)
	{

		filenum++;

		sprintf(weightsfilename, "%sweights_%s.dat", getenv("OUTPUTDIR"), fileext);
		sprintf(vistfilename, "%svis_tot_%s.dat", getenv("OUTPUTDIR"), fileext);
		sprintf(visdfilename, "%svis_diff_%s.dat", getenv("OUTPUTDIR"), fileext);
		sprintf(noisetfilename, "%snoise_tot_%s.dat", getenv("OUTPUTDIR"), fileext);
		sprintf(noisedfilename, "%snoise_diff_%s.dat", getenv("OUTPUTDIR"), fileext);

		/**************************************************/

		w1 = calloc(size_half * Nchan, sizeof(float));

		/* Compute weights contribution */

		if ((fptrv = fopen(weightsfilename, "r")) == NULL)
		{
			// file does not exist - ignore
			fprintf(flog, "Weights file does not exist %s... ignoring...\n", weightsfilename);
		}
		else
		{

			fprintf(flog, "filename: %s\n", weightsfilename);

			if (fread(w1, sizeof(float), size_fin * Nchan, fptrv) != size_fin * Nchan)
				fprintf(flog, "Error: input file is the incorrect size. File: %s\n", weightsfilename);

			fclose(fptrv);

			printf("Adding\n");

			// annoying row-major order
			for (i = 0; i < size_half * Nchan; i++)
				weights[i] += w1[i]; /// norm[ch];
		}

		/**************************************************/

		// Compute other contributions - weighted

		/* Set-up output arrays */

		vist = calloc(size_half * Nchan, sizeof(float complex));
		visd = calloc(size_half * Nchan, sizeof(float complex));
		noiset = calloc(size_half * Nchan, sizeof(float complex));
		noised = calloc(size_half * Nchan, sizeof(float complex));

		/* Compute vist contribution */

		if ((fptrv = fopen(vistfilename, "r")) == NULL)
		{
			// file does not exist - ignore
			fprintf(flog, "File does not exist %s... ignoring...\n", vistfilename);
		}
		else
		{

			fprintf(flog, "filename: %s\n", vistfilename);

			if (fread(vist, sizeof(float complex), size_fin * Nchan, fptrv) != size_fin * Nchan)
				fprintf(flog, "Error: input file is the incorrect size. File: %s\n", vistfilename);

			fclose(fptrv);
		}

		/***********************/

		/* Compute visd contribution */

		if ((fptrv = fopen(visdfilename, "r")) == NULL)
		{
			// file does not exist - ignore
			fprintf(flog, "File does not exist %s... ignoring...\n", visdfilename);
		}
		else
		{

			fprintf(flog, "filename: %s\n", visdfilename);

			if (fread(visd, sizeof(float complex), size_fin * Nchan, fptrv) != size_fin * Nchan)
				fprintf(flog, "Error: input file is the incorrect size. File: %s\n", visdfilename);

			fclose(fptrv);
		}

		/***********************/

		/* Compute noiset contribution */

		if ((fptrv = fopen(noisetfilename, "r")) == NULL)
		{
			// file does not exist - ignore
			fprintf(flog, "File does not exist %s... ignoring...\n", noisetfilename);
		}
		else
		{

			fprintf(flog, "filename: %s\n", noisetfilename);

			if (fread(noiset, sizeof(float complex), size_fin * Nchan, fptrv) != size_fin * Nchan)
				fprintf(flog, "Error: input file is the incorrect size. File: %s\n", noisetfilename);

			fclose(fptrv);
		}

		/***********************/

		/* Compute noised contribution */

		if ((fptrv = fopen(noisedfilename, "r")) == NULL)
		{
			// file does not exist - ignore
			fprintf(flog, "File does not exist %s... ignoring...\n", noisedfilename);
		}
		else
		{

			fprintf(flog, "filename: %s\n", noisedfilename);

			if (fread(noised, sizeof(float complex), size_fin * Nchan, fptrv) != size_fin * Nchan)
				fprintf(flog, "Error: input file is the incorrect size. File: %s\n", noisedfilename);

			fclose(fptrv);
		}

		/***********************/

		/****** Perform the weighted sum *********/

		//    if ((addsub == 0)||((filenum == 1)&&(addsub == 1))){

		/* Subtract the mean first */

		fprintf(flog, "Adding\n");
		for (i = 0; i < size_half * Nchan; i++)
		{
			vistot[i] += vist[i] * w1[i];
			visdiff[i] += visd[i] * w1[i];
			noisetot[i] += noiset[i] * w1[i];
			noisediff[i] += noised[i] * w1[i];
		}

		/*    } else {


				fprintf(flog,"Subtracting\n");

				for (i=0;i<size_half*Nchan;i++){
					vistot[i] -= vist[i]*w1[i];
					visdiff[i] -= visd[i]*w1[i];
					noisetot[i] -= noiset[i]*w1[i];
					noisediff[i] -= noised[i]*w1[i];
				}


			}
		*/
		fprintf(flog, "weights %f %f %f visdiff %f %f\n", w1[15000], (w1[15001]), (w1[15002]), creal(visd[150]), cimag(visd[150]));

		/**************** Free internal arrays ************/

		free(vist);
		free(visd);
		free(noiset);
		free(noised);
		free(w1);

	} // End loop through list of input data to combine

	fclose(fin);

	/***************************************************************************************/
	/******************* Re-weight outputs before writing **********************************/

	fprintf(flog, "Total weights %f vistot %f %f visdiff %f %f\n", weights[150], creal(vistot[150]), cimag(vistot[150]), creal(visdiff[150]), cimag(visdiff[150]));

	for (i = 0; i < size_half * Nchan; i++)
	{
		if (weights[i] != 0.)
		{
			vistot[i] /= weights[i];
			visdiff[i] /= weights[i];
			noisetot[i] /= weights[i];
			noisediff[i] /= weights[i];
		}
		else
		{
			vistot[i] = 0.;
			visdiff[i] = 0.;
			noisetot[i] = 0.;
			noisediff[i] = 0.;
		}
	}

	fprintf(flog, "Total post division weights %f vistot %f %f visdiff %f %f\n", weights[150], creal(vistot[150]), cimag(vistot[150]), creal(visdiff[150]), cimag(visdiff[150]));

	// write_out data for visibility contributions

	sprintf(filename_real2, "%svis_diff_%s.dat", getenv("OUTPUTDIR"), ext);
	sprintf(filename_real, "%svis_tot_%s.dat", getenv("OUTPUTDIR"), ext);

	fptrv = fopen(filename_real2, "w");
	if (fptrv == NULL)
	{
		printf("Error: can't open output file of vis.\n");
	}
	else
	{
		printf("File opened successfully.\n");
		for (i = 0; i < size_half * Nchan; i++)
		{
			//      double temp=creal(out_array[i]);
			fwrite(&(visdiff[i]), sizeof(float complex), 1, fptrv);
		}
	}
	fclose(fptrv);

	fptrv = fopen(filename_real, "w");
	if (fptrv == NULL)
	{
		printf("Error: can't open output file of vis.\n");
	}
	else
	{
		printf("File opened successfully.\n");
		for (i = 0; i < size_half * Nchan; i++)
		{
			//      double temp=creal(out_array[i]);
			fwrite(&(vistot[i]), sizeof(float complex), 1, fptrv);
		}
	}
	fclose(fptrv);

	// write_out data for noise contributions

	sprintf(bvfilename, "%snoise_diff_%s.dat", getenv("OUTPUTDIR"), ext);
	sprintf(bbfilename, "%snoise_tot_%s.dat", getenv("OUTPUTDIR"), ext);

	fptrv = fopen(bvfilename, "w");
	if (fptrv == NULL)
	{
		printf("Error: can't open output file of noise.\n");
	}
	else
	{
		printf("File opened successfully.\n");
		for (i = 0; i < size_half * Nchan; i++)
		{
			//      double temp=creal(out_array[i]);
			fwrite(&(noisediff[i]), sizeof(float complex), 1, fptrv);
		}
	}
	fclose(fptrv);

	fptrv = fopen(bbfilename, "w");
	if (fptrv == NULL)
	{
		printf("Error: can't open output file of noise.\n");
	}
	else
	{
		printf("File opened successfully.\n");
		for (i = 0; i < size_half * Nchan; i++)
		{
			//      double temp=creal(out_array[i]);
			fwrite(&(noisetot[i]), sizeof(float complex), 1, fptrv);
		}
	}
	fclose(fptrv);

	// write_out data for weights contribution

	sprintf(normfilename, "%sweights_%s.dat", getenv("OUTPUTDIR"), ext);

	fptrv = fopen(normfilename, "w");
	if (fptrv == NULL)
	{
		printf("Error: can't open output file of weights.\n");
	}
	else
	{
		printf("File opened successfully.\n");
		for (i = 0; i < size_half * Nchan; i++)
		{
			fwrite(&(weights[i]), sizeof(float), 1, fptrv);
		}
	}
	fclose(fptrv);

	/* write completion to log */
	fprintf(flog, "Completed.\n");

	fclose(flog);

	return 0;
}

// Functions

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
