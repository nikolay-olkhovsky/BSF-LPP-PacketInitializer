/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP starting point Initializer
Module: Problem-Forwards.h (Problem Function Forwards)
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-bsfTypes.h"
#include "Problem-Types.h"
//====================== User Functions ==============================
void packetGenerator_Init(bool* success);
void packetGenerator_Run();

//====================== Problem Functions ===========================
bool		Like(PT_vector_T a1, PT_float_T b1, PT_float_T a1Norm, PT_vector_T a2, PT_float_T b2, PT_float_T a2Norm);
int			RndSign();
PT_float_T	RndValue(PT_float_T rndMax);
void		RndVector(PT_vector_T vector);
void		Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint);
void		Vector_MultiplyByNumber(PT_vector_T x, PT_float_T r, PT_vector_T y);
PT_float_T	Vector_NormSquare(PT_vector_T x);
void		Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z);

#ifdef PP_DATABASE_OUTPUT
std::vector<double> charToDouble(std::vector<char> _In);
std::vector<char> doubleToChar(std::vector<double> _In);
void printLppForm(Problem problem, std::vector<Inequality> inequalities);
#endif //PP_DATABASE_OUTPUT
//====================== Macros ================================
#define PF_MIN(x,y) (x<y?x:y)
#define PF_MAX(x,y) (x>y?x:y)