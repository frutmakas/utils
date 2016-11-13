#ifndef _INTERPOL_H_
#define _INTERPOL_H_


void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

void ratint(double xa[], double ya[], int n, double x, double *y, double *dy);

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);


#endif

