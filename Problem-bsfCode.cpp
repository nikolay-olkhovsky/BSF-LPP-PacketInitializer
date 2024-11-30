/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP Generator
Module: Problem-bsfCode.cpp (Implementation of the Problem)
Prefix: PC
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
// PP_STATE_START
// PP_STATE_FIND_INITIAL_APPROXIMATION
// PP_STATE_FIND_INTERIOR_POINT
// PP_STATE_DETERMINE_DIRECTION
// PP_STATE_MOVING_ALONG_SURFACE
// PP_STATE_LANDING

#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	//???
};

void PC_bsf_Start(bool* success) {
	ini::IniFile config;

	cout << "PC_bsf_Start beginnig!" << endl;
	if (PP_BSF_MAX_MPI_SIZE != 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "PP_BSF_MAX_MPI_SIZE must be equal to 2!" << endl;
		*success = false;
		return;
	}


	config.load(PP_FILE_INI);
	PP_PATH = config["general"]["PP_PATH"].as<string>();
	PP_PROBLEM_NAME = config["general"]["PP_PROBLEM_NAME"].as<string>();
	PP_MTX_PREFIX = config["general"]["PP_MTX_PREFIX"].as<string>();
	PP_MTX_POSTFIX_A = config["general"]["PP_MTX_POSTFIX_A"].as<string>();
	PP_MTX_POSTFIX_B = config["general"]["PP_MTX_POSTFIX_B"].as<string>();
	PP_MTX_POSTFIX_LO = config["general"]["PP_MTX_POSTFIX_LO"].as<string>();
	PP_MTX_POSTFIX_HI = config["general"]["PP_MTX_POSTFIX_HI"].as<string>();
	PP_MTX_POSTFIX_C = config["general"]["PP_MTX_POSTFIX_C"].as<string>();
	PP_MTX_POSTFIX_X0 = config["general"]["PP_MTX_POSTFIX_X0"].as<string>();
	PP_MTX_POSTFIX_SO = config["general"]["PP_MTX_POSTFIX_SO"].as<string>();
	PP_LPP_FILE = config["general"]["PP_LPP_FILE"].as<string>();

	PP_N = config["general"]["PP_N"].as<int>();
	PP_NUM_OF_RND_INEQUALITIES = config["generator"]["PP_NUM_OF_RND_INEQUALITIES"].as<int>();
	PP_RND_SEED = config["generator"]["PP_RND_SEED"].as<int>();
	PP_MTX_N = (2 * PP_N + PP_NUM_OF_RND_INEQUALITIES + 1); // ???
	if(PP_NUM_OF_RND_INEQUALITIES == 0)
		PP_M = (PP_N + 1);
	else
		PP_M = (PP_N + PP_NUM_OF_RND_INEQUALITIES);
	PP_MTX_M = (PP_N + PP_NUM_OF_RND_INEQUALITIES + 1); // ???
	PP_ALPHA = config["generator"]["PP_ALPHA"].as<float>();
	
	PP_THETA = (PP_ALPHA / 2);
	PP_RHO = (PP_THETA / 2);
	PP_A_MAX = (PP_RHO / 2);
	PP_B_MAX = (PP_ALPHA * 50);
	PP_MAX_LIKE = config["generator"]["PP_MAX_LIKE"].as<float>();
	PP_LIKE_FACTOR = config["generator"]["PP_LIKE_FACTOR"].as<float>();
	PP_MIN_SHIFT = (PP_RHO / 3);

	PP_NUMBER_OF_PROBLEMS = config["generator"]["PP_NUMBER_OF_PROBLEMS"].as<unsigned long>();
	PP_OUTPUT_LIMIT = config["generator"]["PP_OUTPUT_LIMIT"].as<unsigned long>();
	PP_SETW = config["generator"]["PP_SETW"].as<unsigned long>();

	PD_index = 0;
	PD_time = clock();
	PD_seed = PD_previous_seed = 0;
	//srand(time(0));

	if(BSF_sv_mpiRank == BSF_sv_mpiMaster)
	{
#ifdef PP_DATABASE_OUTPUT
		storage.sync_schema();
		storage.begin_transaction();
#else
		PD_packetWriter = new CMTXPacketWriter(PP_PATH.c_str(), PP_NUMBER_OF_PROBLEMS);
		PD_packetWriter->clearFolder();
		PD_packetWriter->open();
		PD_lppPacketWriter = new CMTXLppPacketWriter(PP_PATH.c_str(), PP_NUMBER_OF_PROBLEMS);
		PD_lppPacketWriter->open();
#endif // PP_DATABASE_OUTPUT
	}
	PD_failuresType1 = 0;
	PD_failuresType2 = 0;
	cout << "PC_bsf_Start ending!" << endl;
}

void PC_bsf_Init(bool* success) {
	if (BSF_sv_mpiRank == BSF_sv_mpiMaster) {
		if (PP_RND_SEED > 0)
			PD_seed = PP_RND_SEED;
		else
			while (PD_seed == PD_previous_seed) {
				PD_seed = (unsigned)clock() * (BSF_sv_mpiRank + 10);
			}
		//cout << PD_seed << endl;
		srand(PD_seed);
		PD_previous_seed = PD_seed;
		packetGenerator_Init(success);
		packetGenerator_Run();
	}
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = BSF_sv_numOfWorkers;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	//elem->a = PD_A[i];
	//elem->b = &(PD_b[i]);
}

// 0. Pseudo-pojection
void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success // 1 - reduceElem was produced successfully; 0 - otherwise
) {
	/*
	PT_float_T aNormSquare;
	PT_float_T term; // computing b
	bool failure = false;
	bool like = false;
	
	reduceElem->failuresType1 = 0;
	reduceElem->failuresType2 = 0;
	
	if (PP_NUM_OF_RND_INEQUALITIES == 0)
		return;

	do {
		for (int j = 0; j < PD_n; j++) // computing a[*]
			reduceElem->a[j] = PP_THETA + RndSign() * RndValue(PP_THETA);
		term = PP_RHO / PD_sqrt_n + PP_ALPHA / 2 + RndValue((PP_THETA - PP_RHO) / PD_sqrt_n);
		reduceElem->b = 0;
		for (int j = 0; j < PD_n; j++)
			reduceElem->b += reduceElem->a[j] * term;

		failure = false;
		for (int j = 0; j < PD_n; j++)
			if (reduceElem->b / reduceElem->a[j] <= PP_ALPHA) { // Point (0,...,0,PP_ALPHA,0,...,0) is not feasible
				failure = true;
				break;
			}
		if (failure) {
			reduceElem->failuresType2++;
			continue;
		}

		aNormSquare = Vector_NormSquare(reduceElem->a);
		reduceElem->aNorm = sqrt(aNormSquare);

		like = false;
		for (int i = 0; i < PD_m_predef; i++)
			if (like = Like(reduceElem->a, reduceElem->b, reduceElem->aNorm, PD_A[i], PD_b[i], PD_aNorm[i]))
				break;

		if (like) {// Hyperplane is similar to another one
			reduceElem->failuresType1++;
			continue;
		}

		break;
	} while (true);
	*/
}

// 1. CheckIn
void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {
	// not used
}

void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {
	// not used
}

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {
	// not used
}

// 0. Pseudo-pojection
void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	// z = x (+) y
	//z->failuresType1 = x->failuresType1 + y->failuresType1;
	//z->failuresType2 = x->failuresType2 + y->failuresType2;
}

// 1. CheckIn
void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {
	// not used
}

void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {
	// not used
}

void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {
	// not used
}

//0. Start
void PC_bsf_ProcessResults(
	PT_bsf_reduceElem_T* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	bool success; // to call PC_bsf_Init(bool* success)
	if (PD_index < PP_NUMBER_OF_PROBLEMS) {
	//if (PD_index < 2) {
		*nextJob = BD_JOB_RESET;
		return;
	}
	else {
#ifndef PP_DATABASE_OUTPUT
		PD_packetWriter->close();
		PD_lppPacketWriter->close();
#else
		storage.commit();
#endif //PP_DATABASE_OUTPUT
		*exit = true;
		return;
	}
}

// 1. Movement on Polytope  ========================================================
void PC_bsf_ProcessResults_1(
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
}

void PC_bsf_ProcessResults_2(
	PT_bsf_reduceElem_T_2* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
}

void PC_bsf_ProcessResults_3(
	PT_bsf_reduceElem_T_3* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
}

void PC_bsf_JobDispatcher(
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* job,
	bool* exit,
	double t
) {
	// Optional filling. Do not delete!
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
//	cout << "-------------------------------------PC_bsf_ParametersOutput-----------------------------------" << endl;
//	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
//#ifdef PP_BSF_OMP
//#ifdef PP_BSF_NUM_THREADS
//	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
//#else
//	cout << "Number of Threads: " << omp_get_num_procs() << endl;
//#endif // PP_BSF_NUM_THREADS
//#else
//	cout << "OpenMP is turned off!" << endl;
//#endif // PP_BSF_OMP
//	cout << "Dimension n = " << PP_N << endl;
//	cout << "Number of random inequalities: " << PP_NUM_OF_RND_INEQUALITIES << endl;
//	cout << "Length of hypercube edge ALPHA = " << PP_ALPHA << endl;
//	cout << "Radius of large hypersphere RHO = " << PP_THETA << endl;
//	cout << "Radius of small hypersphere THETA = " << PP_RHO << endl;
//	cout << "Maximal acceptable likeness of equations MAX_LIKE = " << PP_LIKE_FACTOR << endl;
//	cout << "Minimal acceptable shift MIN_SHIFT = " << PP_MIN_SHIFT << endl;
//#ifdef PP_MATRIX_OUTPUT
//	cout << "------- Support inequalities -------" << endl;
//	for (int i = 0; i < PD_m; i++) {
//		cout << i << ")";
//		for (int j = 0; j < PD_n; j++)
//			cout << setw(PP_SETW) << PD_A[i][j] << "\t";
//		cout << "\t<=\t" << setw(PP_SETW) << PD_b[i] << endl;
//	}
//	cout << "n = " << PD_n << "\tm = " << PD_m << "\tk = " << PD_k << endl;
//#endif // PP_MATRIX_OUTPUT
//
//	cout << "Objective Function:\t";
//	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(PP_SETW) << PD_c[j];
//	cout << (PP_OUTPUT_LIMIT < PP_N ? "	..." : "") << endl;
}

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	// not used
}

// 0. Start
void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {	// For Job 0
	static int k = PD_m;
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;

	if (PD_k > k) {
		cout << PD_k - 1 << ")\t";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_A[PD_k - 1][j] << "\t";
		cout << (PP_OUTPUT_LIMIT < PD_n ? "	..." : "") << "<=\t" << setw(PP_SETW) << PD_b[PD_k - 1] << endl;
		k = PD_k;
	}
	cout << "Failures 'Similar' = " << PD_failuresType1 << endl;
	cout << "Failures '(0," << PP_ALPHA << ",0)' = " << PD_failuresType2 << endl;
	cout << "-------------------------------------" << endl;
}

// 1. Movement on Polytope
void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	/*cout << "------------------ 1. Movement on Polytope. Iter # " << BSF_sv_iterCounter << " ------------------" << endl;
	cout << "Elapsed time: " << round(elapsedTime) << endl;
	cout << "PD_u:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_u[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
	cout << "PD_direction:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_direction[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
	cout << "Sift Length = " << PD_shiftLength << endl;/**/
};

// 2.
void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
}


void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
}

// 0. Start
void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	if (PD_index % 1000 == 0)
		cout << PD_index << " problems written. Failures = "<< PD_failuresType1 + PD_failuresType2 <<", Time: " << (clock() - PD_time)/CLOCKS_PER_SEC << endl;
//	cout << "=============================================" << endl;
//	cout << "Time: " << t << endl;
//	cout << "Iterations: " << BSF_sv_iterCounter << endl;
//#ifdef PP_MATRIX_OUTPUT
//	cout << "------- Random inequalities -------" << endl;
//	for (int i = PD_m_predef; i < PD_m; i++) {
//		cout << i << ")\t";
//		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_A[i][j] << "\t";
//		cout << (PP_OUTPUT_LIMIT < PD_n ? "	..." : "") << "<=\t" << setw(PP_SETW) << PD_b[i] << endl;
//	}
//	cout << "-----------------------------------" << endl;
//#endif // PP_MATRIX_OUTPUT
//	cout << "Failures 'Similar' = " << PD_failuresType1 << endl;
//	cout << "Failures '(0," << PP_ALPHA << ",0)' = " << PD_failuresType2 << endl;
}

// 1. Movement on Polytope
void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	//ProblemOutput(t);
}

void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// not used
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// not used
}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; };
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; };
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; };
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; };
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; };
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; };
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; };
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; };

//---------------------------------- User functions ----------------------------
inline void packetGenerator_Init(bool* success) {
#ifndef PP_DATABASE_OUTPUT
	PD_currentProblem = new CProblem;
	PD_currentProblem->dimensions = PP_N;
	PD_currentProblem->constraints = (PP_NUM_OF_RND_INEQUALITIES < 1 ? 1 : PP_NUM_OF_RND_INEQUALITIES);
#endif //PP_DATABASE_OUTPUT
	PD_centerObjectF = 0;
	PD_k = 2 * PP_N + 1;

	//srand((unsigned)time(NULL) * (BSF_sv_mpiRank + 10));
	PD_n = PP_N;

	if (PD_n < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "PP_N must be greater than 1!" << endl;
		*success = false;
		return;
	}

	PD_m = PD_n;
	PD_k = PD_n;

#ifndef PP_DATABASE_OUTPUT
	PD_currentProblem->A = new CMatrix(PP_M, PP_N, 0);
	PD_currentProblem->lpp = new CMatrix(PP_M + PP_N, PP_N, 0);
	PD_currentProblem->b = new CArray(PP_M + PP_N);
#endif //PP_DATABASE_OUTPUT
	for (int i = 0; i < PP_N; i++) {
		for (int j = 0; j < PP_N; j++)
		{
			PD_A[i][j] = 0;
#ifndef PP_DATABASE_OUTPUT
			PD_currentProblem->A->setValue(i, j, 0.);
			PD_currentProblem->lpp->setValue(i, j, 0.);
#endif //PP_DATABASE_OUTPUT
		}
		PD_A[i][i] = 1;
#ifndef PP_DATABASE_OUTPUT
		PD_currentProblem->A->setValue(i, i, 1.);
		PD_currentProblem->lpp->setValue(i, i, 1.);
#endif //PP_DATABASE_OUTPUT
//		Vector_Copy(PD_A[i], PD_dataset[PD_index].A[i]);
		PD_b[i] = PP_ALPHA;
//		PD_dataset[PD_index].b[i] = PD_b[i];
#ifndef PP_DATABASE_OUTPUT
		PD_currentProblem->b->setValue(i, PD_b[i]);
#endif //PP_DATABASE_OUTPUT
	}

	if (PP_NUM_OF_RND_INEQUALITIES == 0)
	{
		for (int j = 0; j < PP_N; j++)
		{
			PD_A[PP_N][j] = 1;
#ifndef PP_DATABASE_OUTPUT
			PD_currentProblem->A->setValue(PP_N, j, 1.);
			PD_currentProblem->lpp->setValue(PP_N, j, 1.);
#endif //PP_DATABASE_OUTPUT
		}
//		Vector_Copy(PD_A[PP_N], PD_dataset[PD_index].A[PP_N]);
		PD_b[PP_N] = PP_ALPHA * (PP_N - 1) + PP_ALPHA / 2;
//		PD_dataset[PD_index].b[PP_N] = PD_b[PP_N];
#ifndef PP_DATABASE_OUTPUT
		PD_currentProblem->b->setValue(PP_N, PD_b[PP_N]);
#endif //PP_DATABASE_OUTPUT
		PD_m++; assert(PD_m <= PP_MTX_M);
		PD_k++; assert(PD_k <= PP_MTX_M);
	}

	PD_m_predef = PD_m;
	/*
		for (int i = PD_m_predef; i < PD_m_predef + PP_N; i++) {
			for (int j = 0; j < PP_N; j++)
			{
				//PD_A[i][j] = 0;
				//PD_currentProblem->A->setValue(i, j, 0.);
				PD_currentProblem->lpp->setValue(i, j, 0.);
			}
			//PD_A[i][i - PP_N - 1] = -1;
			//PD_currentProblem->A->setValue(i, i - PP_N - 1, -1.);
			PD_currentProblem->lpp->setValue(i, i - PP_N - 1, -1.);
			//Vector_Copy(PD_A[i], PD_dataset[PD_index].A[i]);
			//PD_b[i] = 0;
			//PD_dataset[PD_index].b[i] = PD_b[i];
			PD_currentProblem->b->setValue(i, 0.);
		}
	*/
	for (int j = 0; j < PP_N; j++)
		PD_center[j] = PP_ALPHA / 2;

#ifndef PP_DATABASE_OUTPUT
	PD_currentProblem->c = new CArray(PP_N);
	PD_currentProblem->hi = new CArray(PP_N);
	PD_currentProblem->lo = new CArray(PP_N);
#endif //PP_DATABASE_OUTPUT
	for (int j = 0; j < PP_N; j++) {
		if (PP_NUM_OF_RND_INEQUALITIES == 0) // Standart objective function
			PD_c[j] = (PT_float_T)(PP_N - j) * PP_RHO;
		else // Random objective function
			PD_c[j] = (PT_float_T)(rand() % PD_n + 1);
//		PD_dataset[PD_index].c[j] = PD_c[j];
		PD_centerObjectF += PD_c[j] * PP_THETA;

#ifndef PP_DATABASE_OUTPUT
		PD_currentProblem->c->setValue(j, PD_c[j]);
		PD_currentProblem->hi->setValue(j, 1.0e+308);
		PD_currentProblem->lo->setValue(j, 0.);
#endif //PP_DATABASE_OUTPUT
	}

	for (int i = 0; i < PD_m; i++)
		PD_aNorm[i] = sqrt(Vector_NormSquare(PD_A[i]));

	PD_sqrt_n = (PT_float_T)sqrt(PD_n);
}

void packetGenerator_Run() {
	PT_vector_T a;
	PT_float_T b;
	PT_float_T aNorm;

	PT_float_T aNormSquare;
	PT_float_T term; // computing b
	bool failure = false;
	bool like = false;

	if (PD_k == PP_NUM_OF_RND_INEQUALITIES + PD_m_predef) {
#ifdef PP_DATABASE_OUTPUT
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster) {
			Problem problem;
			Inequality inequality;
			std::vector<Inequality> inequalities;
			std::vector<double> _buff(PP_N);

			problem.N = PP_N;
			problem.seed = PD_seed;
			problem.high = PP_ALPHA;
			problem.low = 0.;
			for (unsigned i = 0; i < PP_N; i++)
				_buff[i] = PD_c[i];
			problem.c = doubleToChar(_buff);
			problem.id = storage.insert<Problem>(problem);

			inequalities.clear();
			for (unsigned i = 0; i < PP_NUM_OF_RND_INEQUALITIES; i++) {
				inequality.problem_id = problem.id;
				inequality.b = PD_b[PP_N + i];
				_buff.resize(PP_N);
				for (unsigned j = 0; j < PP_N; j++)
					_buff[j] = PD_A[PP_N + i][j];
				inequality.coefficients = doubleToChar(_buff);
				inequalities.push_back(inequality);
			}
			storage.insert_range<Inequality>(inequalities.begin(), inequalities.end());
		}
#else
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			PD_packetWriter->addProblem(*PD_currentProblem->convertToMTX());

		for (int i = PP_M; i < PP_M + PP_N; i++) {
			for (int j = 0; j < PP_N; j++)
				PD_currentProblem->lpp->setValue(i, j, 0.);
			PD_currentProblem->lpp->setValue(i, i - PP_N - 1, -1.);
			PD_currentProblem->b->setValue(i, 0.);
		}

		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			PD_lppPacketWriter->addProblem(*PD_currentProblem);
#endif //PP_DATABASE_OUTPUT
		PD_index++;
		return;
	}

	do {
		for (int j = 0; j < PD_n; j++) // computing a[*]
			a[j] = PP_THETA + RndSign() * RndValue(PP_THETA);
		term = PP_RHO / PD_sqrt_n + PP_ALPHA / 2 + RndValue((PP_THETA - PP_RHO) / PD_sqrt_n);
		b = 0;
		for (int j = 0; j < PD_n; j++)
			b += a[j] * term;

		failure = false;
		for (int j = 0; j < PD_n; j++)
			if (b / a[j] <= PP_ALPHA) { // Point (0,...,0,PP_ALPHA,0,...,0) is not feasible
				failure = true;
				break;
			}
		if (failure) {
			PD_failuresType2++;
			continue;
		}

		aNormSquare = Vector_NormSquare(a);
		aNorm = sqrt(aNormSquare);

		like = false;
		for (int i = 0; i < PD_m_predef; i++)
			if (like = Like(a, b, aNorm, PD_A[i], PD_b[i], PD_aNorm[i]))
				break;

		if (like) {// Hyperplane is similar to another one
			PD_failuresType1++;
			continue;
		}

		Vector_Copy(a, PD_A[PD_k]);
//		Vector_Copy(PD_A[PD_k], PD_dataset[PD_index].A[PD_k]);
#ifndef PP_DATABASE_OUTPUT
		for (unsigned i = 0; i < PP_N; i++)
		{
			PD_currentProblem->A->setValue(PD_k, i, PD_A[PD_k][i]);
			PD_currentProblem->lpp->setValue(PD_k, i, PD_A[PD_k][i]);
		}
#endif //PP_DATABASE_OUTPUT
		PD_b[PD_k] = b;
//		PD_dataset[PD_index].b[PD_k] = PD_b[PD_k];
		PD_aNorm[PD_k] = aNorm;
#ifndef PP_DATABASE_OUTPUT
		PD_currentProblem->b->setValue(PD_k, PD_b[PD_k]);
#endif //PP_DATABASE_OUTPUT
		PD_k++;

		if (PD_k == PP_NUM_OF_RND_INEQUALITIES + PD_m_predef) {
#ifdef PP_DATABASE_OUTPUT
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster) {
				Problem problem;
				Inequality inequality;
				std::vector<Inequality> inequalities;
				std::vector<double> _buff(PP_N);
			
				problem.N = PP_N;
				problem.seed = PD_seed;
				problem.high = PP_ALPHA;
				problem.low = 0.;
				for (unsigned i = 0; i < PP_N; i++)
					_buff[i] = PD_c[i];
				problem.c = doubleToChar(_buff);
				problem.id = storage.insert<Problem>(problem);

				inequalities.clear();
				for (unsigned i = 0; i < PP_NUM_OF_RND_INEQUALITIES; i++) {
					inequality.problem_id = problem.id;
					inequality.b = PD_b[PP_N + i];
					_buff.resize(PP_N);
					for (unsigned j = 0; j < PP_N; j++)
						_buff[j] = PD_A[PP_N + i][j];
					inequality.coefficients = doubleToChar(_buff);
					inequalities.push_back(inequality);
				}
				storage.insert_range<Inequality>(inequalities.begin(), inequalities.end());
			}
#else
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster) {
				PD_packetWriter->addProblem(*PD_currentProblem->convertToMTX());
			}

			for (int i = PP_M; i < PP_M + PP_N; i++) {
				for (int j = 0; j < PP_N; j++)
					PD_currentProblem->lpp->setValue(i, j, 0.);
				PD_currentProblem->lpp->setValue(i, i - PP_N - 1, -1.);
				PD_currentProblem->b->setValue(i, 0.);
			}
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				PD_lppPacketWriter->addProblem(*PD_currentProblem);
#endif //PP_DATABASE_OUTPUT
			PD_index++;
			return;
		}
	} while (true);
}

//---------------------------------- Problem functions -------------------------
inline PT_float_T Vector_NormSquare(PT_vector_T x) {
	PT_float_T s = 0;

	for (int j = 0; j < PD_n; j++)
		s += x[j] * x[j];
	return s;
}

inline void RndVector(PT_vector_T vector) {
	for (int j = 0; j < PD_n; j++)
		vector[j] = RndValue(PP_A_MAX);
}

inline PT_float_T RndValue(PT_float_T rndMax) { // rnd >= 0
	return ((PT_float_T)rand() / ((PT_float_T)RAND_MAX + 1)) * rndMax;
}

inline int RndSign() {
	int res = rand() % 2;
	if (res == 0)
		res = -1;
	return res;
}

inline bool Like(PT_vector_T a1, PT_float_T b1, PT_float_T a1Norm, PT_vector_T a2, PT_float_T b2, PT_float_T a2Norm) {
	PT_float_T like, shift;
	PT_vector_T e1, e2, e1_e2;

	Vector_MultiplyByNumber(a1, 1 / a1Norm, e1);
	Vector_MultiplyByNumber(a2, 1 / a2Norm, e2);

	Vector_Subtraction(e1, e2, e1_e2);
	like = sqrt(Vector_NormSquare(e1_e2));
	shift = fabs(b1 / a1Norm - b2 / a2Norm);
	if (like < PP_LIKE_FACTOR)
		if (shift < PP_MIN_SHIFT)
			return true;
	return false;
}

inline void Vector_MultiplyByNumber(PT_vector_T x, PT_float_T r, PT_vector_T y) {  // y = r*x
	for (int j = 0; j < PD_n; j++) {
		y[j] = x[j] * r;
	}
}

inline void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
	for (int j = 0; j < PD_n; j++) {
		z[j] = x[j] - y[j];
	}
}

inline void Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
	for (int j = 0; j < PD_n; j++) {
		toPoint[j] = fromPoint[j];
	}
}

#ifdef PP_DATABASE_OUTPUT
std::vector<double> charToDouble(std::vector<char> _In) {
	std::vector<double> _Out;
	size_t size = _In.size();
	_Out.resize(size / sizeof(double));
	std::memcpy(_Out.data(), _In.data(), size);
	return _Out;
}

std::vector<char> doubleToChar(std::vector<double> _In) {
	std::vector<char> _Out;
	size_t size = _In.size() * sizeof(double);
	_Out.resize(size);
	std::memcpy(_Out.data(), _In.data(), size);
	return _Out;
}

void printLppForm(Problem problem, std::vector<Inequality> inequalities) {
	int M = 2 * problem.N + inequalities.size();
	int width = 10;
	std::vector<double> _vec;
	std::cout << problem.id << '\t' << problem.N << '\t' << M << std::endl;
	for (int i = 0; i < problem.N; i++) {
		for (int j = 0; j < problem.N; j++)
			if (i == j) std::cout << std::setw(width) << double(1);
			else std::cout << std::setw(width) << double(0);
		std::cout << std::setw(width) << problem.high << std::endl;
	}
	for (int i = 0; i < inequalities.size(); i++) {
		_vec = charToDouble(inequalities[i].coefficients);
		for (int j = 0; j < problem.N; j++)
			std::cout << std::setw(width) << _vec[j];
		std::cout << std::setw(width) << inequalities[i].b << std::endl;
	}
	for (int i = 0; i < problem.N; i++) {
		for (int j = 0; j < problem.N; j++)
			if (i == j) std::cout << std::setw(width) << double(-1);
			else std::cout << std::setw(width) << double(0);
		std::cout << std::setw(width) << problem.low << std::endl;
	}
	_vec = charToDouble(problem.c);
	std::copy(_vec.begin(), _vec.end(), std::ostream_iterator<double>(std::cout, "\t"));
	std::cout << std::endl;
}
#endif // PP_DATABASE_OUTPUT

