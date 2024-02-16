/*==============================================================================
Project: LiFe
Theme: Packet LPP Generator
Module: Problem-bsfParameters.h (LBSF-skeleton Parameters)
Prefix: PP_BSF
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton
==============================================================================*/

//=========================== BSF Skeleton Parameters =========================
#define PP_BSF_MAX_MPI_SIZE 500		// Defines the maximum possible number of MPI processes
#define PP_BSF_PRECISION 6			// Sets the decimal precision to be used to format floating-point values on output operations
//#define PP_BSF_ITER_OUTPUT			// If it is defined then Iteration Output is performed
#define PP_BSF_TRACE_COUNT 3000		// Each PP_BSF_TRACE_COUNT-th iteration to be outputted
#define PP_BSF_MAX_JOB_CASE 0		// Defines the maximum number of activities (jobs) in workflow minus 1

#define BD_JOB_RESET -1				// Special job to reinitialize start data for a task

//--------------------------- OpenMP Parameters ---------------------------
//#define PP_BSF_OMP				// If PP_BSF_OMP is defined then OpenMP is turned on for Map Step
//#define PP_BSF_NUM_THREADS 6	// If PP_BSF_NUM_THREADS is udefined then all accessable threads are used