/*******************************************************************************************
MAP decoder Simulation in C

Note :
-------
This program is used to simulate the Maximum A-Posteriori Probability algorithm in C.
Since this program is used only for simulation, I have not done much of code optimization (as
well as memory optimization). So this program might run a bit slower. I'll try to optimize the
code and memory requirements and update this program as soon as possible.

Parameters :
---------------
Thus MAP algorithm will decode a single convolutional code. The parameters of the convolutional
code can be easily set by modifying the N_STATES, POSSIBLE_INPUTS, NextStates[] and ActualOutputs[][]
in accordance with the code that you have used.

Documentation :
----------------
Detailed comments are provided throughout this Code. There will not be much problems in reading it.
For some good explanation about this algorithm, please visit this page.

http://www.redrival.com/rvc/faq/map.htm

This page is a FAQ on this algorithm that I have prepared. It is better to read that FAQ first before
going thro' this code.
For comments and suggestion drop me line at vc@lucent.com.

Important Notes :
------------------
1. I have used Likehood ratios instead of Log-Likehood ratios currently, because of some numerical
   instability problems. I will correct the code later.
2. The trellis is assumed to be terminated. If this is not the case, the funtion InitializeBetas()
   has to be modified. (A general purpose function can also be provided.. but unnecessary 'cause
   the trellis is normally terminated in most cases).


********************************************************************************************/

#include <stdio.h>
#include <math.h>
#include "intleave.h"

/* Some global Constants */
/* Defines the possible no of inputs for a given state. It will be 2^(input block length).*/

#define POSSIBLE_INPUTS 2
#define N_STATES 4 /* Total no of states in the state machine */
#define BLOCK_LENGTH 18 /* Length of the block to be decoded */
#define NUM_OUTPUTS 2 /* No of outputs in the Convolutional Encoder. */

/* List of Global Variables and Arrays used */
/* ======================================== */
/* This table holds the values of the next states */
/* This table can also be calculated dynamically using the static variables and bit shifts.*/
int NextStates[N_STATES*POSSIBLE_INPUTS] = {0,2,2,0,3,1,1,3};
/* This array holds the actual outputs of the channel encoder */
int ActualOutputs[N_STATES*POSSIBLE_INPUTS][NUM_OUTPUTS] = {{0,0},{1,1},{0,0},{1,1},{0,1},{1,0},{0,1},{1,0}};
/* This array holds tha outputs that were received from the demodulator (Corrupted Outputs) */
//int  NoisyOutputs[NUM_OUTPUTS][BLOCK_LENGTH] = {{0,0},{1,1},{0,1},{1,0},{0,1},{1,1},{0,0},{1,1},{0,1},{1,0},{0,1},{1,1},{0,0},{1,1},{0,1},{1,0},{0,1},{1,1}};
int NormalOutputs[BLOCK_LENGTH][NUM_OUTPUTS] = {{0,0},{1,1},{0,1},{1,0},{0,1},{1,1},{0,0},{1,1},{0,1},{1,0},{0,1},{1,1},{0,0},{1,1},{0,1},{1,0},{0,1},{1,1}};
/* This table will hold the values of alpha */
double Alpha[N_STATES][BLOCK_LENGTH+1];
/* This table will hold the values of beta */
double Beta[N_STATES][BLOCK_LENGTH+1];
/* This array will hold the outputs of the Channel Demodulator */
/* Since they are declared as 'int', they can be configured for soft-decision */
int Outputs[BLOCK_LENGTH][NUM_OUTPUTS];
/* This array will hold the resultant LLR's after the decoding process. */
/* These LLR's will constitute the input to the next stage */
/* Please note that, I am currently using just LR's instead of LLR's */
double LLR[BLOCK_LENGTH];

/* This table holds the apriori probabilities (from the previous decoder or demodulator)        */
/* Remember that this probabilities are Log-Likelihood ratios in the case of iterative decoding */
double aPrioriProbs[BLOCK_LENGTH] = {.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,0.5,0.5};

/* Function Declarations. -- You won't see all the declarations here */
/* Please see the function implementation for details of its execution and interfaces */
int GetNextState(int PresentState,int EncoderInputs);
double GetAprioriProbability(int TimeIndex);

//-----------------------------------------------------------------------------------------//
// GetNextState(..) - Function
//==============================
// This function gets the next state for a given present state..
// currently uses table look-up... can be converted such that the corresponding
// values can be calculated using bit-shifts.
//
// Inputs:
// =======
// PresentState -- This is the present state bit vector of the state machine.
// EncoderInputs -- This is the input bit vector to the encoder state machine
//
// Output:
// =======
// The Next State bit vector of the state machine.
//-----------------------------------------------------------------------------------------//
int GetNextState(int PresentState, int EncoderInputs)
{
	return NextStates[PresentState*POSSIBLE_INPUTS+EncoderInputs];
}

//-----------------------------------------------------------------------------------------//
// GetAprioriProbability(..) - Function
// =====================================
// This function returns the apriori probabilities of the input at the given
// time index. The LLR's are converted into probabilities in this function. The
// Probability of being a '0' is returned. The other prob can be found by subtracting
// it from 1.
//
// Note :: Currently I am using just the Likelihood ratios, rather than log-likelihood
// ratios (LLR's gave some math instability problems..).
//
// Inputs
// -------
// int TimeIndex -- The time index at which these probabilities are needed.
//
// Outputs
// --------
// double -- The apriori probability of the input at the given time index being a '0'
//-----------------------------------------------------------------------------------------//
double GetAprioriProbability(int TimeIndex)
{
	return aPrioriProbs[TimeIndex];

//	return ((double)1.0/(1+exp(aPrioriProbs[TimeIndex])));
}

//-----------------------------------------------------------------------------------------//
// IntializeAlphas() - Function
// ============================
// This function Initializes the values of the alpha array in the time index '0'.
// During this intialization we assume that we start in the state number '0'.
//-----------------------------------------------------------------------------------------//
void InitializeAlphas()
{
	int j=0;
	Alpha[0][0] = 1;
	for(j=1 ; j<N_STATES ; j++)
		Alpha[j][0] = 0;
}

//-----------------------------------------------------------------------------------------//
// CalculateEuclideanDistance(..) - Function
// ==========================================
// This function calculates the EuclideanDistance between the DecoderInputs and the
// actual outputs.
//
// Inputs :
// --------
// TimeIndex -- The Decoder Output Index.
// StateNumber, InputOfThatState -- Uniquely identify the actual output that must be
// compared with the decoder output.
//-----------------------------------------------------------------------------------------//
int CalculateEuclideanDistance(int TimeIndex, int StateNumber, int InputOfThatState)
{
	int Sum = 0,i,Value;
	for(i=0;i<NUM_OUTPUTS;i++)
	{
		Value = (NormalOutputs[TimeIndex][i] - ActualOutputs[StateNumber*POSSIBLE_INPUTS+InputOfThatState][i]);
		Sum = Sum + Value*Value;
	}
	return Sum;
}

//-----------------------------------------------------------------------------------------//
// CalcGamma(..) -- Function
// ===========================
// Calculates the value of gamma for the given state transition s(i)->s(i+1).
// This function does the following things :
// -- Get the Apriori probability of the message bit at the given time instant.
// -- Calculate the Euclidean Channel Metric between the inputs and outputs of the channel.
// -- Using these, finds the value of Gamma for the given state transition.
//
// Inputs :
// ---------
// TimeIndex -- The value of 'i'. (Remember, Gamma(s(i)->s(i+1) is to be found).
// State -- the Intial State. (s(i)).
// InputOfThatState -- the Input of that state that will provide the final state s(i+1).
// Note : NextTimeIndex is assumed to be TimeIndex+1.
//-----------------------------------------------------------------------------------------//
double CalcGamma(int TimeIndex, int State, int InputOfThatState)
{
	double AprioriProbability, EuclideanMetric, ChannelMetric;
	// The code that follows works only for the case in which the input block size is 1
	// I'll have to change the code later to include generic block sizes.
	AprioriProbability = GetAprioriProbability(TimeIndex);
	// The AprioriProbability returned by the above function is for '0' only.
	// You need to change it for an input of '1'
	if (InputOfThatState == 1) AprioriProbability = 1 - AprioriProbability;

	// Calculate the channel metric (Euclidean Distance between noisy and actual sequence).
	EuclideanMetric = CalculateEuclideanDistance(TimeIndex, State , InputOfThatState);
	// Channel Metric is actually the exponent of the Euclidean Distance.
	ChannelMetric = exp(-EuclideanMetric);

	// Calculate Gamma then..
	return (AprioriProbability*ChannelMetric);
}

/* This global variable is used to remember the search for previous connected states */
/* See the function GetPreviousConnectedState(..) */
int PrevIndexLeft = 0;
//-----------------------------------------------------------------------------------------//
// void ResetPreviousConnectedStateSearch(..) - Function.
// ======================================================
// Resets the previous search for the connected state in the previous time instant.
// I have to do all these, 'cause there is no proper way to find the previous states
// of a given state.
//-----------------------------------------------------------------------------------------//
void ResetPreviousConnectedStateSearch(void)
{
	PrevIndexLeft = -1;
}

//-----------------------------------------------------------------------------------------//
// GetPreviousConnectedState(..) - Function
// =========================================
// This funtion returns the previous connected state to a given state. The search
// for the previous state resumes from the state next to PrevIndexLeft, so that
// you gotta be careful, if you are starting a new search. (Do a reset always before
// starting a new search.)
//
// Inputs :
// ---------
// CurrentState -- The state whose previous states are being searched.
//
// Outputs :
// ---------
// 1. returns the previous state number
// 2. Sets the 'InputMessageBits' to the corresponding value.
//-----------------------------------------------------------------------------------------//
int GetPreviousConnectedState(int CurrentState, int *InputMessageBits)
{
	int i;
	// Search for the previous connected state in the trellis.
	// Start from the state that is next to the one that was left
	// in the previous search.
	for(i=PrevIndexLeft+1; i<(N_STATES*POSSIBLE_INPUTS); i++)
	{
		if (NextStates[i] == CurrentState)
		{
			PrevIndexLeft = i;
			*InputMessageBits = i%POSSIBLE_INPUTS;
			return i/POSSIBLE_INPUTS;
		}
	}
	// Else if all the states are exhausted, no more previous connected states remain..
	// Send an End of State Message
	PrevIndexLeft = 0;
	return -1;
}

//-----------------------------------------------------------------------------------------//
// UpdateAlpha(..) -- Function
// ===============================
// This function updates the Alpha array for the given time index.
//
// Input :
// -------
// This function updates the Alpha array according to the following equation..
//   		Alpha(j,i) = Sum { Alpha(j',i-1)Gamma(s(i-1)->s(i) }
// where the summation is over all the state transitions between the states in the previous
// time instant and the states in the current time instant, which are connected in the Trellis.
//
// The things that must be done are ::
// ------------------------------------
// For all the states
// -- Reset the Previous State Search Engine.
// -- Initialize the Alpha for the corresponding state (and time instant) to zero.
// -- While the search for the previous state does not give an End of State Mesage
// -- -- Calculate Gamma
// -- -- Update Alpha using the given information above.
// -- Exit While
// Exit For.
// That's all
//----------------------------------------------------------------------------------------//
void UpdateAlpha(int CurTimeIndex)
{
	int i;
	int PreviousState, InputMessageBits;
	int EuclideanMetric;
	double AprioriProbability, ChannelMetric, Gamma;

	for(i=0;i<N_STATES;i++)
	{
		Alpha[i][CurTimeIndex] = 0;
		// Reset the previous state search.. (to start afresh for this new state).
		ResetPreviousConnectedStateSearch();
		while ((PreviousState = GetPreviousConnectedState(i, &InputMessageBits))!=-1)
		{
			// Calculate the value of Gamma for this state transition.
			Gamma = CalcGamma(CurTimeIndex-1, PreviousState, InputMessageBits);
			// Update Alpha with the new value...
			Alpha[i][CurTimeIndex] += Alpha[PreviousState][CurTimeIndex-1]*Gamma;
		}
	}
}

//----------------------------------------------------------------------------------------//
// ForwardRecursion(..) -- Function.
// ==================================
// The forward recursion is used to calculate the Alphas.
// This calls the UpdateAlpha() routine to set the values of all the alphas  ..
// See UpdateAlpha() function for more details.
//----------------------------------------------------------------------------------------//
void ForwardRecursion(void)
{
	int i;
	// First of all, intialize the Alpha array .. for the first time instant only.
	InitializeAlphas();
	// Loop for all the time instants.. till the end of the Trellis is reached.
	for(i=1;i<=BLOCK_LENGTH;i++)
	{
		// Update the Alphas for the current time instant.
		UpdateAlpha(i);
	}
}

//----------------------------------------------------------------------------------------//
// IntializeBetas(..) -- Function.
// ===============================
// This function Initializes the values of the Beta array in the time index 'L'.
// ('L' here is referred as BLOCK_LENGTH.)
// Note :: During this intialization we assume that the trellis is terminated.
//----------------------------------------------------------------------------------------//
void InitializeBetas(void)
{
	int j;
	// Assuming the Trellis is terminated....
	Beta[0][BLOCK_LENGTH] = 1;
	// Initialize the others to zeros..
	for(j=1;j<N_STATES;j++)
	{
		Beta[j][BLOCK_LENGTH] = 0;
	}
}


//----------------------------------------------------------------------------------------//
// UpdateBeta(..) -- Function..
// =============================
// This function updates the Beta array for the given time index.
//
// Input :
// -------
// This function updates the Beta array according to the following equation..
//   		Beta(j,i) = Sum { Beta(j',i+1)Gamma(s(i)->s(i+1) }
// where the summation is over all the state transitions between the states in the previous
// time instant and the states in the current time instant, which are connected in the Trellis.
//
// The things that must be done are ::
// ------------------------------------
// For all the states
// -- Reset the Previous State Search Engine.
// -- Initialize the Beta value for the corresponding state (and time instant) to zero.
// -- While the search for the previous state does not give an End of State Mesage
// -- -- Calculate Gamma
// -- -- Update Beta using the given information above.
// -- Exit While
// Exit For.
// That's all...
//----------------------------------------------------------------------------------------//
void UpdateBeta(int CurTimeIndex)
{
	// Here 'i' is used as a state index and j is used as an input for the index.
	int i,j;
	// To hold the Next State Variable ..
	int NextState, EuclideanMetric;
	double AprioriProbability, Gamma, ChannelMetric;

	// Loop over all the states in the Trellis at the "CurTimeIndex" instant..
	for(i=0;i<N_STATES;i++)
	{
		Beta[i][CurTimeIndex] = 0;
		// Loop over all possible inputs..
		// Here, we are a bit lucky, as there is no searching  business here.
		// Everything is straight-forward memory reference.
		for(j=0;j<POSSIBLE_INPUTS;j++)
		{
			// Get the next state...
			NextState = GetNextState(i,j);
			// Get the value of Gamma for this state transition...
			Gamma = CalcGamma(CurTimeIndex,i,j);
			// Update the corresponding Beta array..
			Beta[i][CurTimeIndex] += (Beta[NextState][CurTimeIndex+1]*Gamma);
		}
	}
}

//----------------------------------------------------------------------------------------//
// BackwardRecursion(..) -- Function
// ==================================
// The Backward Recursion is used to calculate the Betas. This is done by calling the
// UpdateBeta() function for every time instant beginning from the last
// Please see UpdateBeta() functio header for details.
//----------------------------------------------------------------------------------------//
void BackwardRecursion(void)
{
	int i;
	// Initialize the Beta array .. for the last time instant only.
	InitializeBetas();
	// Update the Betas array for all the time instants .. starting from the last.
	for(i=BLOCK_LENGTH-1;i>=0;i--)
	{
		UpdateBeta(i);
	}
}

//-----------------------------------------------------------------------------------------//
// ComputeLLRs(..) -- Function
// ============================
// This is the last stage of MAP decoding. It computes the Likelihood ratios.
// See the MAP FAQ for details of how this is done.
//-----------------------------------------------------------------------------------------//
void ComputeLLRs(void)
{
	int i,j;
	// Accumulators for Prob[0] and Prob[1].
	double Sum0, Sum1;
	// NextState Variables..
	int NextState0, NextState1;
	// For all the time instants...
	for(i=0;i<BLOCK_LENGTH;i++)
	{
		// Reset the accumulators every time instant..
		Sum0 = 0; Sum1 = 0;
		// For all the states in the given time instant...
		for(j=0;j<N_STATES;j++)
		{
			// Get the two possible next states for inputs 0 and 1.
			NextState0 = GetNextState(j,0);
			NextState1 = GetNextState(j,1);

			// Calculate the sums using the formulae for LR's. (See MAP FAQ..)
			Sum0 += Alpha[j][i]*CalcGamma(i,j,0)*Beta[NextState0][i+1];
			Sum1 += Alpha[j][i]*CalcGamma(i,j,1)*Beta[NextState1][i+1];
		}
		// Calculate LLR.
		LLR[i] = Sum1/Sum0;
	}
}

//-----------------------------------------------------------------------------------------//
// MAPDecode(..) -- Function
// ==========================
// Doesnt do much except calling ForwardRecursion(), BackwardRecursion() and ComputeLLRs().
//-----------------------------------------------------------------------------------------//
void MAPDecode(void)
{
	ForwardRecursion();
	BackwardRecursion();
	ComputeLLRs();
}

//-----------------------------------------------------------------------------------------//
// This is the "main" main().
//-----------------------------------------------------------------------------------------//
#define NUM_ITERS 5 // Number of iterations for the MAP algorithm..
void main(void)
{
	int j;
	// Execute the MAP Algorithm
	MAPDecode();
	printf("Results of the MAP Algorithm ::\n");
	for(j=0;j<BLOCK_LENGTH;j++)
	{
		printf("LR[%d] = %lf (Likelihood - %d)\n",j,LLR[j],(LLR[j]>1)? 1 : 0);
	}
}


