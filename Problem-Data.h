/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP Generator
Module: Problem-Data.h (Problem Data)
Prefix: PD
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-Types.h"		// Problem Parameters 
using namespace std;

static PT_unsigned_T PP_NUMBER_OF_PROBLEMS;
static PT_unsigned_T PP_OUTPUT_LIMIT;
static PT_unsigned_T PP_SETW;

static PT_float_T PP_ALPHA;								// Length of hypercube edge
static PT_float_T PP_THETA;								// Radius of large hypersphere
static PT_float_T PP_RHO;								// Radius of small hypersphere
static PT_float_T PP_A_MAX;								// Maximal random value for A
static PT_float_T PP_B_MAX;								// Maximal random value for b
static PT_float_T PP_MAX_LIKE;							// Maximal acceptable likeness of equations (must be less then 0.7)
static PT_float_T PP_LIKE_FACTOR;
static PT_float_T PP_MIN_SHIFT;							// Minimal acceptable shift

//========================== Problem variables ====================================
static int PD_n;					// Current dimension
static int PD_m;					// Current number of inequalities
static int PD_m_predef;				// Number of predefined inequalities
static PT_float_T PD_sqrt_n;		// Square root of n
static PT_float_T PD_centerObjectF;	// Value of object function in the center of hypercube
static unsigned PD_failuresType1 = 0;
static unsigned PD_failuresType2 = 0;
static unsigned PD_failuresType3 = 0;
static int PD_k;	// Index of current random inequality
static int PD_seed;
static int PD_previous_seed;
static PT_unsigned_T PD_index = 0;			// Index of current LPP in dataset
static double PD_time;						// Overall time of execution

//========================== Problem data structures ==============================
static PT_matrix_T PD_A;
static PT_column_T PD_b;
static PT_vector_T PD_c;
static PT_vector_T PD_center;		// Center of hypercube
static PT_column_T PD_aNorm;
vector<PT_MTXrow_T> PD_MTXdataset(0);

//========================== Files ==============================
#ifdef PP_DATABASE_OUTPUT
static auto storage = sqlite_orm::make_storage("C:/HS/dataset100000_4.sqlite3",
	sqlite_orm::make_table("problems",
		sqlite_orm::make_column("id", &Problem::id, sqlite_orm::primary_key()),
		sqlite_orm::make_column("N", &Problem::N),
		sqlite_orm::make_column("seed", &Problem::seed),
		sqlite_orm::make_column("high", &Problem::high),
		sqlite_orm::make_column("low", &Problem::low),
		sqlite_orm::make_column("c", &Problem::c)
	),
	sqlite_orm::make_table("inequalities",
		sqlite_orm::make_column("id", &Inequality::id, sqlite_orm::primary_key()),
		sqlite_orm::make_column("coefficients", &Inequality::coefficients),
		sqlite_orm::make_column("b", &Inequality::b),
		sqlite_orm::make_column("problem_id", &Inequality::problem_id),
		sqlite_orm::foreign_key(&Inequality::problem_id).references(&Problem::id)
	),
    sqlite_orm::make_table("surface_points",
        sqlite_orm::make_column("id", &SurfacePoint::id, sqlite_orm::primary_key()),
        sqlite_orm::make_column("coefficients", &SurfacePoint::coefficients),
        sqlite_orm::make_column("problem_id", &SurfacePoint::problem_id),
        sqlite_orm::foreign_key(&SurfacePoint::problem_id).references(&Problem::id)
    ),
    sqlite_orm::make_table("precedents",
        sqlite_orm::make_column("id", &Precedent::id, sqlite_orm::primary_key().autoincrement()),
        sqlite_orm::make_column("coefficients", &Precedent::coefficients),
        sqlite_orm::make_column("face", &Precedent::face),
        sqlite_orm::make_column("d", &Precedent::d),
        sqlite_orm::make_column("shift", &Precedent::shift),
        sqlite_orm::make_column("face_count", &Precedent::face_count),
        sqlite_orm::make_column("face_numbers", &Precedent::face_numbers),
        sqlite_orm::make_column("problem_id", &Precedent::problem_id),
        sqlite_orm::foreign_key(&Precedent::problem_id).references(&Problem::id)
    ),
    sqlite_orm::make_table("images",
        sqlite_orm::make_column("id", &Image::id, sqlite_orm::primary_key().autoincrement()),
        sqlite_orm::make_column("density", &Image::density),
        sqlite_orm::make_column("rank", &Image::rank),
        sqlite_orm::make_column("answer_vector", &Image::answer_vector),
        sqlite_orm::make_column("cosine_vector", &Image::cosine_vector),
        sqlite_orm::make_column("num_of_points", &Image::num_of_points),
        sqlite_orm::make_column("data", &Image::data),
        sqlite_orm::make_column("precedent_id", &Image::precedent_id),
        sqlite_orm::make_column("field_points", &Image::field_points),
        sqlite_orm::foreign_key(&Image::precedent_id).references(&Precedent::id)
    )
);
#else
CProblem* PD_currentProblem;
PT_dataset_T* PD_dataset = new PT_dataset_T[PP_PROBLEMS_LIMIT];		// Many LPPs in one dataset
CMTXPacketWriter*		PD_packetWriter;
CMTXLppPacketWriter*	PD_lppPacketWriter;
#endif // PP_DATABASE_OUTPUT

//static string PD_fileName;

// nor - number of rows; noc - number of columns; non - number of non-zero values
/*static FILE* PD_stream_A;
static string PD_MTX_File_A; /* File of matrix A (only non-zero elements):
------------ begin of file -------------
nor noc non // nor=m (number of inequalities); noc=n (dimension); non - number of non-zero values
i_1 j_1 A_{{i_1}{j_1}} // i_1 - index of row; j_1 - index of column
...
i_k j_k A_{{i_k}{j_k}}
------------ end of file----------------*/
/*static FILE* PD_stream_b;
static string PD_MTX_File_b; /* File of column b:
------------ begin of file -------------
nor noc // nor=m (number of inequalities); noc=1
b_1
...
b_{nor}
------------ end of file----------------*/
/*static FILE* PD_stream_lo;
static string PD_MTX_File_lo; /* File of variable lower bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
lo_1
...
lo_{nor}
------------ end of file----------------*/
/*static FILE* PD_stream_hi;
static string PD_MTX_File_hi; /* File of variable higher bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
lo_1
...
lo_{nor}
------------ end of file----------------*/
/*static FILE* PD_stream_c;
static string PD_MTX_File_c; /* File of variable higher bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
c_1
...
c_{nor}
------------ end of file----------------*/
/*static FILE* PD_stream_x0;
static string PD_MTX_File_x0;/* File of starting point in the following format:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
x_1
...
x_{nor}
------------ end of file----------------*/
/*static string PD_MTX_File_so;/* Solution file in the following format:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
x_1
...
x_{nor}
------------ end of file----------------*/
//static FILE* PD_stream_tr;
//static string PD_MTX_File_tr;