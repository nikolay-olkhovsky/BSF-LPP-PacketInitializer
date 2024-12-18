/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP Generator
Module: Problem-Types.h (LBSF Types)
Prefix: PT
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/			
#pragma once
#include "Problem-Include.h"		// Problem "Include" Files
#include "Problem-Parameters.h"		// Problem Parameters 
//=========================== Problem Types =========================
typedef float PT_float_T;
typedef unsigned long long PT_unsigned_T;
typedef PT_float_T	PT_vector_T[PP_MAX_N];
typedef PT_float_T	PT_matrix_T[PP_MAX_M][PP_MAX_N];
typedef PT_float_T	PT_column_T[PP_MAX_M];
typedef struct {
	PT_matrix_T A;
	PT_column_T b;
	PT_vector_T c;
} PT_dataset_T;
typedef struct {
	int row;
	int col;
	PT_float_T val;
} PT_MTXrow_T;

typedef PT_float_T	PT_MTXvector_T[PP_MAX_MTX_N];
typedef PT_float_T	PT_MTXmatrix_T[PP_MAX_MTX_M][PP_MAX_MTX_N];
typedef PT_float_T	PT_MTXcolumn_T[PP_MAX_MTX_M];

#ifdef PP_DATABASE_OUTPUT
struct Problem {
    unsigned id;
    int N;
    int seed;
    double high;
    double low;
    std::vector<char> c;
};

struct Inequality {
    unsigned id;
    std::vector<char> coefficients;
    double b;
    int problem_id;
};

struct SurfacePoint {
    unsigned id;
    std::vector<char> coefficients;
    int problem_id;
};

struct Precedent {
    unsigned id;
    int problem_id;
    std::vector<char> coefficients;
    int face;
    std::vector<char> d;
    double shift;
    std::vector<char> face_count;
    std::vector<char> face_numbers;
};

struct Image {
    unsigned id;
    int precedent_id;
    double density;
    double rank;
    std::vector<char> answer_vector;
    std::vector<char> cosine_vector;
    int num_of_points;
    std::vector<char> data;
    std::optional<std::vector<char>> field_points;
};
#endif // PP_DATABASE_OUTPUT
