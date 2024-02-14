/*==============================================================================
Project: LiFe
Theme: Packet LPP Generator
Module: Problem-bsfTypes.h (Predefined Problem-depended LBSF Types)
Prefix: PT_bsf
Author: Nikolay A. Olkhovsky
This source code is a part of BSF Skeleton
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 
//=========================== BSF Types =========================
struct PT_bsf_parameter_T {		// Type of Parameter for workers (current approximation)
};

struct PT_bsf_mapElem_T {		// Type of map-list elements
};

struct PT_bsf_reduceElem_T {	// Type of reduce-list elements for Job 0 (default)	
	PT_vector_T a;
	PT_float_T b;
	PT_float_T aNorm;
	unsigned failuresType1;
	unsigned failuresType2;
	unsigned failuresType3;
};

struct PT_bsf_reduceElem_T_1 {	// Type of reduce-list elements for Job 1
	// Not used
};

struct PT_bsf_reduceElem_T_2 {	// Type of reduce-list elements for Job 2
	// Not used
};

struct PT_bsf_reduceElem_T_3 {	// Type of reduce-list elements for Job 3
	// Not used
};