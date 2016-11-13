#ifndef _RAYLEIGH_H_
#define _RAYLEIGH_H_

extern NSPDRandGausState rayleighDGausStatePtr;

double *RayleighFading(double doppler_freq, double sample_freq, int nbSamples, int *output_size);

#endif