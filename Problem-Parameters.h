/*==============================================================================
Project: LiFe
Theme: Packet LPP Generator
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
//=========================== Problem Parameters =========================
#define PP_MAX_RND_INEQUALITIES	40
#define PP_MAX_N				20
#define PP_MAX_MTX_N			81
#define PP_MAX_M				81
#define PP_MAX_MTX_M			61
/*
#define PP_N 3												// Space dimension
#define PP_MTX_N (2*PP_N + PP_NUM_OF_RND_INEQUALITIES + 1)	// Number of columns (variables) in MTX format
#define PP_NUM_OF_RND_INEQUALITIES (2 * PP_N)				// Number of random inequalities
//#define PP_NUM_OF_RND_INEQUALITIES 0
#define PP_M (2*PP_N + PP_NUM_OF_RND_INEQUALITIES + 1)		// Total number of inequalities
#define PP_MTX_M (PP_N + PP_NUM_OF_RND_INEQUALITIES + 1)	// Number of inequalities in MTX matrix
#define PP_ALPHA 200										// Length of hypercube edge
#define PP_THETA (PP_ALPHA/2)								// Radius of large hypersphere
#define PP_RHO (PP_THETA/2)									// Radius of small hypersphere
#define	PP_A_MAX (PP_ALPHA*5)								// Maximal random value for A
#define	PP_B_MAX (PP_ALPHA*50)								// Maximal random value for b
#define	PP_MAX_LIKE 0.8										// Maximal acceptable likeness of equations (must be less then 0.7)
#define	PP_MIN_SHIFT (2*PP_RHO)								// Minimal acceptable shift
*/
static int PP_N;												// Space dimension
static int PP_MTX_N;											// Number of columns (variables) in MTX format
static int PP_NUM_OF_RND_INEQUALITIES;							// Number of random inequalities
static int PP_RND_SEED;
static int PP_M;												// Total number of inequalities
static int PP_MTX_M;											// Number of inequalities in MTX matrix

//-------------------------- Input/Outpoot Parameters ---------------------------
#define PP_PROBLEMS_LIMIT 100000
//#define PP_FILE_OUTPUT
#define PP_FILE_INI "config.ini"
static std::string PP_PATH;
static std::string PP_LPP_FILE;

//------------------------- Matrix format ----------------
static std::string PP_PROBLEM_NAME;
static std::string PP_MTX_PREFIX;
static std::string PP_MTX_POSTFIX_A;
static std::string PP_MTX_POSTFIX_B;
static std::string PP_MTX_POSTFIX_LO;
static std::string PP_MTX_POSTFIX_HI;
static std::string PP_MTX_POSTFIX_C;
static std::string PP_MTX_POSTFIX_X0;
static std::string PP_MTX_POSTFIX_SO;
