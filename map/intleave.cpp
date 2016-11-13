#include <stdio.h>
#include <math.h>

// The interleaver dimension. The actual size of the interleaver is the square of this ..
#define INTERLEAVER_DIMENSION 64

void Interleave(int *InputBitStream, int *OutputBitStream)
{
	int i,j;
	for(i=0;i<INTERLEAVER_DIMENSION;i++)
	{
		for(j=0;j<INTERLEAVER_DIMENSION;j++)
		{
			OutputBitStream[j*INTERLEAVER_DIMENSION+i] = InputBitStream[i*INTERLEAVER_DIMENSION+j];
		}
	}
}

void DeInterleave(int *InputBitStream, int *OutputBitStream)
{
	int i,j;
	for(i=0;i<INTERLEAVER_DIMENSION;i++)
	{
		for(j=0;j<INTERLEAVER_DIMENSION;j++)
		{
			OutputBitStream[j*INTERLEAVER_DIMENSION+i] = InputBitStream[i*INTERLEAVER_DIMENSION+j];
		}
	}
}

void InterleaveDouble(double *InputBitStream, double *OutputBitStream)
{
	int i,j;
	for(i=0;i<INTERLEAVER_DIMENSION;i++)
	{
		for(j=0;j<INTERLEAVER_DIMENSION;j++)
		{
			OutputBitStream[j*INTERLEAVER_DIMENSION+i] = InputBitStream[i*INTERLEAVER_DIMENSION+j];
		}
	}
}

void DeInterleaveDouble(double *InputBitStream, double *OutputBitStream)
{
	int i,j;
	for(i=0;i<INTERLEAVER_DIMENSION;i++)
	{
		for(j=0;j<INTERLEAVER_DIMENSION;j++)
		{
			OutputBitStream[j*INTERLEAVER_DIMENSION+i] = InputBitStream[i*INTERLEAVER_DIMENSION+j];
		}
	}
}

/*
This is just a test routine..
void main(void)
{
	int Input[16],Output[16],i;
	for(i=0;i<16;i++)
	{
		Input[i] = i;
	}
	Interleave(Input,Output);
	for(i=0;i<16;i++)
	{
		printf("\nInput[i] - %d, Output[i] - %d",Input[i],Output[i]);
	}
}
*/
