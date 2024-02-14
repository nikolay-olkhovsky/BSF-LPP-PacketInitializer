/*==============================================================================
Project: LiFe
Theme: Packet LPP Generator
Module: Problem-bsfCode.cpp (Implementation of the Problem)
Prefix: PC
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton
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
	PP_MTX_N = (2 * PP_N + PP_NUM_OF_RND_INEQUALITIES + 1);
	PP_M = (2 * PP_N + PP_NUM_OF_RND_INEQUALITIES + 1);
	PP_MTX_M = (PP_N + PP_NUM_OF_RND_INEQUALITIES + 1);
	PP_ALPHA = config["generator"]["PP_ALPHA"].as<float>();
	
	PP_THETA = (PP_ALPHA / 2);
	PP_RHO = (PP_THETA / 2);
	PP_A_MAX = (PP_ALPHA * 5);
	PP_B_MAX = (PP_ALPHA * 50);
	PP_MAX_LIKE = config["generator"]["PP_MAX_LIKE"].as<float>();
	PP_MIN_SHIFT = (2 * PP_RHO);

	PP_NUMBER_OF_PROBLEMS = config["generator"]["PP_NUMBER_OF_PROBLEMS"].as<unsigned long>();
	PP_OUTPUT_LIMIT = config["generator"]["PP_OUTPUT_LIMIT"].as<unsigned long>();
	PP_SETW = config["generator"]["PP_SETW"].as<unsigned long>();

	PD_index = 0;
	srand(time(0));

	*success = OpenDataFiles();
	if (*success == false)
		return;
	cout << "End of PC_bsf_Start!" << endl;
}

void PC_bsf_Init(bool* success) {
	PD_k = 2 * PP_N + 1;

	srand((unsigned)time(NULL) * (BSF_sv_mpiRank + 10));

	for (int j = 0; j < PP_N; j++)
		PD_center[j] = PP_THETA;

	for (int j = 0; j < PP_N; j++) {
		PD_c[j] = (PT_float_T)(PP_N - j) * PP_RHO;
		PD_dataset[PD_index].c[j] = PD_c[j];
		PD_centerObjectF += PD_c[j] * PP_THETA;
	}

	for (int i = 0; i < PP_N; i++) {
		for (int j = 0; j < PP_N; j++)
			PD_A[i][j] = 0;
		PD_A[i][i] = 1;
		Vector_Copy(PD_A[i], PD_dataset[PD_index].A[i]);
		PD_b[i] = PP_ALPHA;
		PD_dataset[PD_index].b[i] = PD_b[i];
	}

	for (int j = 0; j < PP_N; j++)
		PD_A[PP_N][j] = 1;
	Vector_Copy(PD_A[PP_N], PD_dataset[PD_index].A[PP_N]);
	PD_b[PP_N] = PP_ALPHA * (PP_N - 1) + PP_ALPHA / 2;
	PD_dataset[PD_index].b[PP_N] = PD_b[PP_N];

	for (int i = PP_N + 1; i < 2 * PP_N + 1; i++) {
		for (int j = 0; j < PP_N; j++)
			PD_A[i][j] = 0;
		PD_A[i][i - PP_N - 1] = -1;
		Vector_Copy(PD_A[i], PD_dataset[PD_index].A[i]);
		PD_b[i] = 0;
		PD_dataset[PD_index].b[i] = PD_b[i];
	}

	for (int i = 0; i < 2 * PP_N + 1; i++)
		PD_aNorm[i] = sqrt(Vector_NormSquare(PD_A[i]));
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
	PT_float_T distToCenter;
	PT_float_T aNormSquare;
	PT_vector_T centerProjection;
	PT_float_T aDotProductCenter;
	PT_float_T centerProjectionObjectiveF;
	
	reduceElem->failuresType1 = 0;
	reduceElem->failuresType2 = 0;
	reduceElem->failuresType3 = 0;
	
	do {
		RndFloatVector(reduceElem->a);
		if (reduceElem->a[0] == 0) reduceElem->a[0] += (PT_float_T)0.1;
		reduceElem->b = RndFloatValue(PP_B_MAX);

		aNormSquare = Vector_NormSquare(reduceElem->a);
		reduceElem->aNorm = sqrt(aNormSquare);
		aDotProductCenter = Vector_DotProduct(PD_center, reduceElem->a);
		distToCenter = fabs(aDotProductCenter - reduceElem->b) / reduceElem->aNorm;
		if (distToCenter > PP_THETA || distToCenter < PP_RHO) {
			reduceElem->failuresType1++;
			continue;
		}

		ProjectionOnHiperplane(PD_center, reduceElem->a, aNormSquare, aDotProductCenter, reduceElem->b, centerProjection);
		centerProjectionObjectiveF = Vector_DotProduct(PD_c, centerProjection);
		if (PD_centerObjectF > centerProjectionObjectiveF) {
			reduceElem->failuresType2++;
			continue;
		}

		if (!PointIn(PD_center, reduceElem->a, reduceElem->b)) {
			for (int j = 0; j < PP_N; j++)
				reduceElem->a[j] = -reduceElem->a[j];
			reduceElem->b = -reduceElem->b;
		}

		bool like = false;
		for (int i = 0; i < PP_N + 1; i++)
			if (like = Like(reduceElem->a, reduceElem->b, reduceElem->aNorm, PD_A[i], PD_b[i], PD_aNorm[i]))
				break;

		if (like) {
			reduceElem->failuresType3++;
			continue;
		}

		break;
	} while (true);
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
	z->failuresType1 = x->failuresType1 + y->failuresType1;
	z->failuresType2 = x->failuresType2 + y->failuresType2;
	z->failuresType3 = x->failuresType3 + y->failuresType3;
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
		bool like;
		struct extendedReduceElem_T {	// Extended element type of reduce list
			PT_bsf_reduceElem_T elem;	// Element of reduce list
			int reduceCounter;			// Reduce Counter
		};
		extendedReduceElem_T* extendedReduceElem;

		if (PD_k == PP_NUM_OF_RND_INEQUALITIES + 2 * PP_N + 1) {
			PD_k = 2 * PP_N + 1;
			PC_bsf_Init(&success);
			PD_index++;
		}

		PD_failuresType1 += reduceResult->failuresType1;
		PD_failuresType2 += reduceResult->failuresType2;
		PD_failuresType3 += reduceResult->failuresType3;

		extendedReduceElem = (extendedReduceElem_T*)reduceResult;

		/* debug *//*cout << "----------------\n";
		for (int w = 0; w < BSF_sv_numOfWorkers; w++) {
			cout << "w=" << w << ")\t";
			for (int j = 0; j < PP_N; j++) cout << setw(PP_SETW) << extendedReduceElem[w].elem.a[j] << "\t";
			cout << "<=\t" << setw(PP_SETW) << extendedReduceElem[w].elem.b << endl;
		}/* end debug */

		for (int w = 0; w < BSF_sv_numOfWorkers; w++) {

			like = false;

			for (int i = 2 * PP_N + 1; i < PD_k; i++) {
				if (like = Like(extendedReduceElem[w].elem.a, extendedReduceElem[w].elem.b, extendedReduceElem[w].elem.aNorm, PD_A[i], PD_b[i], PD_aNorm[i]))
					break;
			}

			if (like) {
				PD_failuresType3++;
				continue;
			}

			Vector_Copy(extendedReduceElem[w].elem.a, PD_A[PD_k]);
			Vector_Copy(PD_A[PD_k], PD_dataset[PD_index].A[PD_k]);
			PD_b[PD_k] = extendedReduceElem[w].elem.b;
			PD_dataset[PD_index].b[PD_k] = PD_b[PD_k];
			PD_aNorm[PD_k] = extendedReduceElem[w].elem.aNorm;
			PD_k++;

			if (PD_k == PP_NUM_OF_RND_INEQUALITIES + 2 * PP_N + 1) {
				PD_k = 2 * PP_N + 1;
				PC_bsf_Init(&success);
				PD_index++;
			}
		}
	}
	else {
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
	cout << "=================================================== Problem parameters ====================================================" << endl;
	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
#else
	cout << "Number of Threads: " << omp_get_num_procs() << endl;
#endif // PP_BSF_NUM_THREADS
#else
	cout << "OpenMP is turned off!" << endl;
#endif // PP_BSF_OMP
	cout << "Dimension: n = " << PP_N << endl;
	cout << "Number of Random Inegualities: " << PP_NUM_OF_RND_INEQUALITIES << endl;
	cout << "Length of hypercube edge ALPHA = " << PP_ALPHA << endl;
	cout << "Radius of large hypersphere RHO = " << PP_THETA << endl;
	cout << "Radius of small hypersphere THETA = " << PP_RHO << endl;
	cout << "Maximal acceptable likeness of equations MAX_LIKE = " << PP_MAX_LIKE << endl;
	cout << "Minimal acceptable shift MIN_SHIFT = " << PP_MIN_SHIFT << endl;
#ifdef PP_MATRIX_OUTPUT
	cout << "------- Support inequalities -------" << endl;
	for (int i = 0; i < 2 * PP_N + 1; i++) {
		cout << i << ")";
		for (int j = 0; j < PP_N; j++)
			cout << setw(PP_SETW) << PD_A[i][j];
		cout << "\t<=" << setw(PP_SETW) << PD_b[i] << endl;
	}
#endif // PP_MATRIX_OUTPUT

	cout << "Objective Function:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(PP_SETW) << PD_c[j];
	cout << (PP_OUTPUT_LIMIT < PP_N ? "	..." : "") << endl;
}

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	// not used
}

// 0. Start
void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {
	static int k = 2 * PP_N + 1;
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;

	if (PD_k > k) {
		cout << PD_k - 1 << ")\t";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(PP_SETW) << PD_A[PD_k - 1][j] << "\t";
		cout << (PP_OUTPUT_LIMIT < PP_N ? "	..." : "") << "<=\t" << setw(PP_SETW) << PD_b[PD_k - 1] << endl;
		k = PD_k;
	}
	cout << "Failures 'Not between' = " << PD_failuresType1 << endl;
	cout << "Failures 'Obtuse angle to objective' = " << PD_failuresType2 << endl;
	cout << "Failures 'Similar' = " << PD_failuresType3 << endl;
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
	cout << "=============================================" << endl;
	cout << "Time: " << t << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;
#ifdef PP_MATRIX_OUTPUT
	cout << "------- Random inequalities -------" << endl;
	for (int i = 2 * PP_N + 1; i < 2 * PP_N + 1 + PP_NUM_OF_RND_INEQUALITIES; i++) {
		cout << i << ")\t";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(PP_SETW) << PD_A[i][j] << "\t";
		cout << (PP_OUTPUT_LIMIT < PP_N ? "	..." : "") << "<=\t" << setw(PP_SETW) << PD_b[i] << endl;
	}
	cout << "-----------------------------------" << endl;
#endif // PP_MATRIX_OUTPUT
	cout << "Failures 'Not between' = " << PD_failuresType1 << endl;
	cout << "Failures 'Obtuse angle to objective' = " << PD_failuresType2 << endl;
	cout << "Failures 'Similar' = " << PD_failuresType3 << endl;
#ifdef PP_FILE_OUTPUT
	PD_fileName.assign(PP_PATH);
	PD_fileName += PP_LPP_FILE;
	const char* fileName = PD_fileName.c_str();
	FILE* stream;
	cout << "-----------------------------------" << endl;
	stream = fopen(fileName, "w");
	if (stream == NULL) {
		cout << "Failure of opening file " << fileName << "!\n";
		return;
	}
	fprintf(stream, "%lld\n", PP_NUMBER_OF_PROBLEMS);
	for (int index = 0; index < PP_NUMBER_OF_PROBLEMS; index++) {
		fprintf(stream, "%d\t%d\t%d\n", index + 1, PP_M, PP_N);

		for (int i = 0; i < PP_M; i++) {
			for (int j = 0; j < PP_N; j++)
				fprintf(stream, "%.2f\t", PD_dataset[index].A[i][j]);
			fprintf(stream, "%.2f\n", PD_dataset[index].b[i]);
		}
		for (int j = 0; j < PP_N; j++)
			fprintf(stream, "%.2f\t", PD_dataset[index].c[j]);
		fprintf(stream, "\n");
	}
	fclose(stream);
	cout << "LPP is saved into file '" << fileName << "'." << endl;
	cout << "-----------------------------------" << endl;

	// --------- Output to Matrix A file ---------------- //
	cout << "-----------------------------------" << endl;
	fprintf(PD_stream_A, "%lld\n", PP_NUMBER_OF_PROBLEMS);
	for (int index = 0; index < PP_NUMBER_OF_PROBLEMS; index++) {
		ConvertToMTX(index);
		int size = PD_MTXdataset.size();
		fprintf(PD_stream_A, "%d\t%d\t%d\t%d\n", index + 1, PP_MTX_M, PP_MTX_N, size);
		for (int i = 0; i < size; i++)
			fprintf(PD_stream_A, "%d\t%d\t%.1f\n", PD_MTXdataset[i].row, PD_MTXdataset[i].col, PD_MTXdataset[i].val);
	}
	fclose(PD_stream_A);
	cout << "Matrix A is saved into file '" << fileName << "'." << endl;
	cout << "-----------------------------------" << endl;

	// --------- Output to Column b file ---------------- //
	cout << "-----------------------------------" << endl;
	fprintf(PD_stream_b, "%lld\n", PP_NUMBER_OF_PROBLEMS);
	for (int index = 0; index < PP_NUMBER_OF_PROBLEMS; index++) {
		fprintf(PD_stream_b, "%d\t%d\t%d\n", index + 1, PP_MTX_M, 1);
		for (int i = 0; i < PP_M; i++)
			if (i < PP_N + 1 || i > PP_N * 2)
				fprintf(PD_stream_b, "%.1f\n", PD_dataset[index].b[i]);
	}
	fclose(PD_stream_b);
	cout << "Column b is saved into file '" << fileName << "'." << endl;
	cout << "-----------------------------------" << endl;

	// --------- Output to Vector � file ---------------- //
	cout << "-----------------------------------" << endl;
	fprintf(PD_stream_c, "%lld\n", PP_NUMBER_OF_PROBLEMS);
	for (int index = 0; index < PP_NUMBER_OF_PROBLEMS; index++) {
		fprintf(PD_stream_c, "%d\t%d\t%d\n", index + 1, PP_MTX_N, 1);
		for (int i = 0; i < PP_MTX_N; i++)
			if (i < PP_N)
				fprintf(PD_stream_c, "%.1f\n", -PD_dataset[index].c[i]);
			else
				fprintf(PD_stream_c, "%.1f\n", 0.);
	}
	fclose(PD_stream_c);
	cout << "Vector c is saved into file '" << fileName << "'." << endl;
	cout << "-----------------------------------" << endl;

	// --------- Output to HI file ---------------- //
	cout << "-----------------------------------" << endl;
	fprintf(PD_stream_hi, "%lld\n", PP_NUMBER_OF_PROBLEMS);
	for (int index = 0; index < PP_NUMBER_OF_PROBLEMS; index++) {
		fprintf(PD_stream_hi, "%d\t%d\t%d\n", index + 1, PP_MTX_N, 1);
		for (int i = 0; i < PP_MTX_N; i++)
			fprintf(PD_stream_hi, "%e\n", 1e+308);
	}
	fclose(PD_stream_hi);
	cout << "Values HI are saved into file '" << fileName << "'." << endl;
	cout << "-----------------------------------" << endl;

	// --------- Output to LO file ---------------- //
	cout << "-----------------------------------" << endl;
	fprintf(PD_stream_lo, "%lld\n", PP_NUMBER_OF_PROBLEMS);
	for (int index = 0; index < PP_NUMBER_OF_PROBLEMS; index++) {
		fprintf(PD_stream_lo, "%d\t%d\t%d\n", index + 1, PP_MTX_N, 1);
		for (int i = 0; i < PP_MTX_N; i++)
			fprintf(PD_stream_lo, "%.1f\n", 0.);
	}
	fclose(PD_stream_lo);
	cout << "Values LO are saved into file '" << fileName << "'." << endl;
	cout << "-----------------------------------" << endl;

	//// --------- Output to x0 file ---------------- //
	//PD_fileName.assign(PP_PATH);
	//PD_fileName += PP_MTX_PREFIX;
	//PD_fileName += PP_PROBLEM_NAME;
	//PD_fileName += PP_MTX_POSTFIX_X0;
	//fileName = PD_fileName.c_str();

	//cout << "-----------------------------------" << endl;
	//FILE* stream_x0 = fopen(fileName, "w");
	//if (stream_x0 == NULL) {
	//	cout << "Failure of opening file " << fileName << "!\n";
	//	return;
	//}
	//fprintf(stream_x0, "%d\n", PP_NUMBER_OF_PROBLEMS);
	//for (int index = 0; index < PP_NUMBER_OF_PROBLEMS; index++) {
	//	fprintf(stream_x0, "%d\t%d\n", PP_MTX_N, 1);
	//	for (int i = 0; i < PP_MTX_N; i++)
	//		fprintf(stream_x0, "%.1f\n", 0.);
	//}
	//fclose(stream_x0);
	//cout << "Values LO are saved into file '" << fileName << "'." << endl;
	//cout << "-----------------------------------" << endl;

#endif // PP_FILE_OUTPUT
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

//---------------------------------- Problem functions -------------------------
inline bool OpenDataFiles() {
	const char* mtxFile;
	int packageSize = 0;

	//--------------- Opening A ------------------
	PD_MTX_File_A = PP_PATH;
	PD_MTX_File_A += PP_MTX_PREFIX;
	PD_MTX_File_A += PP_PROBLEM_NAME;
	PD_MTX_File_A += PP_MTX_POSTFIX_A;
	mtxFile = PD_MTX_File_A.c_str();
	PD_stream_A = fopen(mtxFile, "w");

	if (PD_stream_A == NULL) {
		//
		cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	//--------------- Opening b ------------------
	PD_MTX_File_b = PP_PATH;
	PD_MTX_File_b += PP_MTX_PREFIX;
	PD_MTX_File_b += PP_PROBLEM_NAME;
	PD_MTX_File_b += PP_MTX_POSTFIX_B;
	mtxFile = PD_MTX_File_b.c_str();
	PD_stream_b = fopen(mtxFile, "w");

	if (PD_stream_b == NULL) {
		//
		cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	//--------------- Opening lo ------------------
	PD_MTX_File_lo = PP_PATH;
	PD_MTX_File_lo += PP_MTX_PREFIX;
	PD_MTX_File_lo += PP_PROBLEM_NAME;
	PD_MTX_File_lo += PP_MTX_POSTFIX_LO;
	mtxFile = PD_MTX_File_lo.c_str();
	PD_stream_lo = fopen(mtxFile, "w");

	if (PD_stream_lo == NULL) {
		//
		cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	//--------------- Opening c ------------------
	PD_MTX_File_c = PP_PATH;
	PD_MTX_File_c += PP_MTX_PREFIX;
	PD_MTX_File_c += PP_PROBLEM_NAME;
	PD_MTX_File_c += PP_MTX_POSTFIX_C;
	mtxFile = PD_MTX_File_c.c_str();
	PD_stream_c = fopen(mtxFile, "w");

	if (PD_stream_c == NULL) {
		//
		cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	//--------------- Opening hi ------------------
	PD_MTX_File_hi = PP_PATH;
	PD_MTX_File_hi += PP_MTX_PREFIX;
	PD_MTX_File_hi += PP_PROBLEM_NAME;
	PD_MTX_File_hi += PP_MTX_POSTFIX_HI;
	mtxFile = PD_MTX_File_hi.c_str();
	PD_stream_hi = fopen(mtxFile, "w");

	if (PD_stream_hi == NULL) {
		//
		cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	return true;
}

inline bool CloseDataFiles() {
	fclose(PD_stream_A);
	fclose(PD_stream_b);
	fclose(PD_stream_c);
	fclose(PD_stream_hi);
	fclose(PD_stream_lo);
	if (PD_stream_x0) fclose(PD_stream_x0);
	if (PD_stream_tr) fclose(PD_stream_tr);
	return true;
}
/*
inline bool SavePoint(PT_vector_T x, const char* filename, double elapsedTime) {
	FILE* stream;

	stream = fopen(filename, "w");
	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << filename << "'.\n";
		return false;
	}

	fprintf(stream, "%c Elapsed time: %.0f\n%d %d\n", '%', elapsedTime, PD_n, 1);

	for (int j = 0; j < PD_n; j++)
		fprintf(stream, "%.14f\n", x[j]);

	fclose(stream);
	return true;
}*/

inline PT_float_T Vector_DotProduct(PT_vector_T x, PT_vector_T y) {
	PT_float_T s = 0;
	for (int j = 0; j < PP_N; j++)
		s += x[j] * y[j];
	return s;
}

inline PT_float_T Vector_NormSquare(PT_vector_T x) {
	PT_float_T s = 0;

	for (int j = 0; j < PP_N; j++)
		s += x[j] * x[j];
	return s;
}

inline void RndFloatVector(PT_vector_T vector) {
	for (int i = 0; i < PP_N; i++)
		vector[i] = RndFloatValue(PP_A_MAX);
}

inline PT_float_T RndFloatValue(PT_float_T rndMax) {
	return RndSign() * (((PT_float_T)rand() / ((PT_float_T)RAND_MAX + 1)) * rndMax);
}

inline int RndSign() {
	int res = rand() % 2;
	if (res == 0)
		res = -1;
	return res;
}

inline bool PointIn(PT_vector_T x, PT_vector_T a, PT_float_T b) { // If the point belonges to the Halfspace <a,x> <= b
	if (Vector_DotProduct(a, x) > b)
		return false;
	else
		return true;
}

inline bool Like(PT_vector_T a1, PT_float_T b1, PT_float_T a1Norm, PT_vector_T a2, PT_float_T b2, PT_float_T a2Norm) {
	PT_float_T like, shift;
	PT_vector_T e1, e2, e1_e2;

	Vector_MultiplyByNumber(a1, 1 / a1Norm, e1);
	Vector_MultiplyByNumber(a2, 1 / a2Norm, e2);

	Vector_Subtraction(e1, e2, e1_e2);
	like = sqrt(Vector_NormSquare(e1_e2));
	shift = fabs(b1 / a1Norm - b2 / a2Norm);
	if (like < PP_MAX_LIKE)
		if (shift < PP_MIN_SHIFT)
			return true;
	return false;
}

inline void ProjectionOnHiperplane(PT_vector_T x, PT_vector_T a, PT_float_T aNormSquare, PT_float_T aDotProductx, PT_float_T b, PT_vector_T projection) {
	PT_float_T fac = (aDotProductx - b) / aNormSquare;
	for (int j = 0; j < PP_N; j++)
		projection[j] = x[j] - fac * a[j];
}

inline void Vector_MultiplyByNumber(PT_vector_T x, PT_float_T r, PT_vector_T y) {  // y = r*x
	for (int j = 0; j < PP_N; j++) {
		y[j] = x[j] * r;
	}
}

inline void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
	for (int j = 0; j < PP_N; j++) {
		z[j] = x[j] - y[j];
	}
}

inline void Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
	for (int j = 0; j < PP_N; j++) {
		z[j] = x[j] + y[j];
	}
}

inline void Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
	for (int j = 0; j < PP_N; j++) {
		toPoint[j] = fromPoint[j];
	}
}

// --------- Conversion dataset to MTX format ----------- //
inline void ConvertToMTX(int index) {
	int row = 0;
	int MTXrow = 1;
	int freeVals = 1;
	PD_MTXdataset.clear();
	while (row < PP_N + 1) {
		for (int j = 0; j < PP_N; j++) {
			PT_float_T value = PD_dataset[index].A[row][j];
			if (value != 0.0) {
				PD_MTXdataset.push_back({ MTXrow, j + 1, value });
			}
		}
		PD_MTXdataset.push_back({ MTXrow, PP_N + freeVals, 1. });
		freeVals++;
		row++;
		MTXrow++;
	}
	row += PP_N;
	while (row < PP_M) {
		for (int j = 0; j < PP_N; j++) {
			PT_float_T value = PD_dataset[index].A[row][j];
			if (value != 0.0) {
				PD_MTXdataset.push_back({ MTXrow, j + 1, value });
			}
		}
		PD_MTXdataset.push_back({ MTXrow, PP_N + freeVals, 1. });
		freeVals++;
		row++;
		MTXrow++;
	}
}
