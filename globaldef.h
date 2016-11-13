/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/09/08 12:09:01 $
 * $Revision: 1.1.1.1.4.33 $
 * $Id: globaldef.h,v 1.1.1.1.4.33 2006/09/08 12:09:01 syed Exp $
 ********************************************************************/
 
#ifndef _GLOBAL_DEF_H_
#define _GLOBAL_DEF_H_


#define INDEX_OUT_OF_RANGE					   1
#define INVALID_RANGE						   2
#define DECOMPOSE_ERROR						   3
#define EMPTY_MATRIX						   4
#define DIVISION_BY_ZERO					   5
#define INVALID_DIMENSION					   6
#define INVALID_ARGUMENT					   7
#define INVALID_OPERATION					   8	

#define EMPTY_VECTOR						3000
#define INVALID_VECTOR						3120
#define INVALID_VECTOR_DIMENSION			3121
#define INCOMPATIBLE_VECTOR_DIMENSION		3322

#define INCOMPATIBLE_MATRIX_DIMENSION		3310
#define INVALID_MATRIX						3311
#define INVALID_MATRIX_DIMENSION			3312
#define MATRIX_INDEX_OUT_OF_RANGE			3313
#define MATRIX_NOT_SQUARE					3314

#define INVALID_OFDM_DATA_SIZE_TO_DECODE	3401
#define INVALID_FFT_SIZE					3402


#define SINGULAR_MATRIX_ERROR				4440
#define SINGULAR_MATRIX_1_ERROR				4441
#define SINGULAR_MATRIX_2_ERROR				4442

#define MEMORY_ALLOCATION_ERROR				5000
#define OUT_OF_MEMORY						5001

#define INVALID_UVECTOR						6000
#define INCONSISTENT_UVECTOR				6001
#define INCONSISTENT_MATRIX_SIZE_IN_UVECTOR 6002
#define EMPTY_UVECTOR						6003

#define BAD_FILE                            7000
#define FILE_NOT_EXIST                      7001
#define FILE_ERROR                          7002

#define TINY								1e-20
#define LARGE								1e20

#define BAD_PARAMETER                       8000

#ifndef WINDOWS
#define WINDOWS
#endif

#if !defined(WINDOWS) && !defined(__unix__) && !defined(DEVCPP)
#define WIN_DOS
#endif

#ifdef DEVCPP
#define __TIMESTAMP__ "N/A in DevCpp"
#endif



//#define LDPC_LLR_INCLUDED

#define NO_FLAT_FADING

//#define USE_INTERLEAVER

//#define NO_LDPC
//#define VERBOSE_DECODING
#endif
