/* cublic spline interpolation. Ripped from Numerical Recipes, complete with
   crap unreadable code, but fixed for normal C array indexing conventions
   RBW: April, 2007 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "cspline.h"

/***********************
 cubic spline- initialiser function
x:   array size n of x values
y:   array size n of y values
yp1: 1st deriv of function at start. Use >= 1.0e30 to signal "natural" spline with zero 2nd deiv at ends 
ypn: 1st deriv of function at end
y2:  resultant 2nd derivs
*************************/
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]) {
  int i,k;
  float p,qn,sig,un,*u;

  u=calloc(n,sizeof(float));
  if (u==NULL) {
    fprintf(stderr,"spline: no malloc for array size %d\n",n);
    exit(-1);
  }

  if (yp1 > 0.99e30) {
    y2[0]=u[0]=0.0; 
  } else { 
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }

  for (i=1;i<n-1;i++) { 
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+2.0;
    y2[i]= (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  } 
  if (ypn > 0.99e30) {
    qn=un=0.0;	/* 2nd derivs are zero at endpoints if derivs are set above 1e30 */
  } else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0); 
  for (k=n-2;k>=0;k--) y2[k] = y2[k]*y2[k+1]+u[k];
  free(u);
}


/*************************
cubic spline (use after calling spline once for a given set of x,y)
xa: array size n of x values
ya: array size n of y values
x:  x value for desired interpolated value
y:  result
**************************/
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y) {
  int klo,khi,k;
  float h,b,a;

  klo = 0;
  khi = n-1;
  while (khi-klo > 1) {
    k=(khi+klo)/2;
    if (xa[k] > x) khi=k;
    else klo=k;
    //     printf("klo %d %d khi\n",klo,khi);
    //  printf("xa %g %g %g %g\n",xa[0],xa[1],xa[2],xa[3]);
  }
  h=xa[khi]-xa[klo];
  if (h == 0.) printf("cspline xa[khi] %g xa[klo] %g\n",xa[khi],xa[klo]);
  assert(h!=0.0);
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


/**********************************
 initialise for a 2D cublic spline
 x1a: array of size m of x values for columns
 x2a: array of size n of x values for rows
 ya:  array of size m pointers to arrays of size n of y values [col][row]
 y2a: array of size m pointers to arrays of size n of y 2nd derivs (result)
 **********************************/
void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a) {
  int j;

  for(j=0; j<m; j++) {
    spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);
  }
}

/**********************************
 2D cubic spline. Basically does a 1D spline on all the rows, then another
 1d spline on the resuling column of interpolates. MUST call splie2 first.
Args:
  x1a: array size m of x values for columns
  x2a: array size n of x values for rows
***********************************/
void splin2(float x1a[], float x2a[], float **ya, float **y2a, int m, int n, float x1, float x2, float *y) {
  int j;
  float *ytmp,*yytmp;

  ytmp  = calloc(m,sizeof(float));
  yytmp = calloc(m,sizeof(float));

  for (j=0;j<m;j++) splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);   /* 1.0e30 signal to use "natural" spline with 2nd deriv zero at ends */
  splint(x1a,yytmp,ytmp,m,x1,y);

  free(yytmp);
  free(ytmp);
}
