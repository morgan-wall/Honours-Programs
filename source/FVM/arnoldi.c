/**
 * @filename arnoldi.c
 * @date 12-11-2014
 * @author Morgan Wall
 *
 * @brief 
 */

#include "matrix.h"
#include <stdlib.h>
#include <math.h>

#include "mex.h"
#include "blas.h"
#include "lapack.h"

#define Q_ROWS m
#define H_ROWS (restart_value + 1)
#define g_ROWS (restart_value + 1)

/* Interface */

typedef enum {ROW, COLUMN} orientation_t;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]);

void arnoldi(mxArray* A, mxArray* L, mxArray* U, double* b, 
	long m, long n, unsigned int restart_value, 
	unsigned int error_tol, mxArray* Q, mxArray* H, mxArray* g);

void copy_array(int n, double* array, double* copy);

void retrieve_array_from_matrix(int rows, int columns, double* matrix, 
	double* copy, int index, orientation_t orientation);

/* Implementation */

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	const unsigned short int INPUT_ARGUMENT_COUNT = 6;
	const unsigned short int OUTPUT_ARGUMENT_COUNT = 3;

	// validate input arguments
	if (nrhs != INPUT_ARGUMENT_COUNT) {
		mexErrMsgTxt("Eight input parameters are required.");
	} else if (nlhs != OUTPUT_ARGUMENT_COUNT) {
		mexErrMsgTxt("Three output parameters are required.");
	}

	// initialise storage parameters
	mxArray *A, *L, *U;
	double* b;
	long m, n;
	unsigned int restart_value, error_tol;

	// retrieve and parse input parameters
	A = (mxArray*)prhs[0];
	L = (mxArray*)prhs[1];
	U = (mxArray*)prhs[2];
	b = mxGetPr(prhs[3]);
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	restart_value = mxGetScalar(prhs[4]);
	error_tol = mxGetScalar(prhs[5]);

	// initialise arnoldi parameters
	mxArray* Q = mxCreateDoubleMatrix((mwSize)m, (mwSize)restart_value, mxREAL);
	mxArray* H = mxCreateDoubleMatrix((mwSize)(restart_value + 1), (mwSize)restart_value, mxREAL);
	mxArray* g = mxCreateDoubleMatrix((mwSize)(restart_value + 1), (mwSize)1, mxREAL);

	arnoldi(A, L, U, b, m, n, restart_value, error_tol, Q, H, g);

	plhs[0] = Q;
	plhs[1] = H;
	plhs[2] = g;
}

void arnoldi(mxArray* A, mxArray* L, mxArray* U, double* b, 
	long m, long n, unsigned int restart_value, 
	unsigned int error_tol, mxArray* Q, mxArray* H, mxArray* g) {

	// initialise mex solver parameters
	int const MEX_OUTPUT_PARMETERS = 1;
	int const MEX_INPUT_PARAMETERS = 2;

	mxArray *lhs[MEX_OUTPUT_PARMETERS], *rhs[MEX_INPUT_PARAMETERS];

	mxArray* tempMxArray = mxCreateDoubleMatrix(m, 1, mxREAL);

	double* q_data = mxGetPr(Q);
	double* h_data = mxGetPr(H);
	double* g_data = mxGetPr(g);

	// initialise solver parameters
	long double_size = sizeof(double);

	double H_diag_temp;
	double g_temp;
	double* temp_basis_vector = mxCalloc(m, double_size);

	// initialise Given's rotations variables
	double* c = mxCalloc(restart_value, double_size);
	double* s = mxCalloc(restart_value, double_size);
	
	// generate initial basis vector
	double* q = mxCalloc(m, double_size);
	double b_norm = dnrm2(&m, b, &double_size);
	for (int i = 0; i < m; i++) {
		q_data[i] = b[i] / b_norm; 
	}

	// generate each basis vector in the Krylov subspace
	for (unsigned int i = 0; i < restart_value; i++) {

		// determine the next Kyrlov subspace basis vector
		retrieve_array_from_matrix(m, restart_value, q_data, 
			temp_basis_vector, i, COLUMN);
		mxSetPr(tempMxArray, temp_basis_vector);

		rhs[0] = L;
		rhs[1] = tempMxArray;

		mexCallMATLAB(MEX_OUTPUT_PARMETERS, lhs, MEX_INPUT_PARAMETERS, rhs, 
			"mldivide");

		rhs[0] = U;
		rhs[1] = lhs[0];

		mexCallMATLAB(MEX_OUTPUT_PARMETERS, lhs, MEX_INPUT_PARAMETERS, rhs, 
			"mldivide");

		rhs[0] = A;
		rhs[1] = lhs[0];

		mexCallMATLAB(MEX_OUTPUT_PARMETERS, lhs, MEX_INPUT_PARAMETERS, rhs, 
			"mtimes");

		q = mxGetPr(lhs[0]);

		// orthogonalise the basis vector
		for (unsigned int j = 0; j <= i; j++) {

			retrieve_array_from_matrix(m, restart_value, q_data, 
				temp_basis_vector, j, COLUMN);
			
			// h_data[i * H_ROWS + j] = 
			// 	ddot(&m, temp_basis_vector, &double_size, q, &double_size);
			h_data[i * H_ROWS + j] = 0;
			for (int k = 0; k < m; k++) {
				h_data[i * H_ROWS + j] = h_data[i * H_ROWS + j] + temp_basis_vector[k] * q[k];
			}

			retrieve_array_from_matrix(Q_ROWS, restart_value, q_data, 
				temp_basis_vector, j, COLUMN);

			// dscal(&m, &h_data[i * H_ROWS + j], temp_basis_vector, &double_size);
			for (int k = 0; k < m; k++) {
				temp_basis_vector[k] = temp_basis_vector[k] * h_data[i * H_ROWS + j];
			}

			for (int k = 0; k < m; k++) {
				q[k] = q[k] - temp_basis_vector[k];
			}
		}

		// // normalise the basis vector
		// h_data[i * H_ROWS + (i+1)] = dnrm2(&m, q, &double_size);
		// for (int j = 0; j < m; j++) {
		// 	q_data[(i + 1) * Q_ROWS + j] = q[i] / h_data[i * H_ROWS + (i+1)];
		// }

		// apply Given's rotations
		for (int j = 0; j <= i; j++) {

			// construct rotations for new vector is Hessenberg matrix
			if (j == i) {
				s[i] = h_data[i * H_ROWS + (i+1)] 
					/ sqrt(pow(h_data[i * H_ROWS + i], 2) + pow(h_data[i * H_ROWS + (i+1)], 2));
				c[i] = h_data[i * H_ROWS + i] 
					/ sqrt(pow(h_data[i * H_ROWS + i], 2) + pow(h_data[i * H_ROWS + (i+1)], 2));
			}

			H_diag_temp = c[j] * h_data[i * H_ROWS + j] 
				+ s[j] * h_data[i * H_ROWS + (j+1)];
			h_data[i * H_ROWS + (j+1)] = - s[j] * h_data[i * H_ROWS + j]
				+ c[j] * h_data[i * H_ROWS + (j+1)];
			h_data[i * H_ROWS + j] = H_diag_temp;
		}

		// // hack (avoid round-off error)
		// h_data[i * H_ROWS + (i+1)] = 0.0;

		g_temp = g_data[i];
		g_data[i] = c[i] * g_temp;
		g_data[i+1] = - s[i] * g_temp;

		if (abs(g_data[i+1]) < error_tol || i == restart_value) {
			break;
		}

		// if (abs((*g)[i+1]) < error_tol || i == restart_value) {

		// 	// // reallocate basis vectors based on the number of Arnoldi iterations
		// 	// for (int j = 0; j < m; j++) {
		// 	// 	(*Q)[j] = realloc((*Q)[j], sizeof(***Q) * i);
		// 	// }

		// 	// // reallocate Hessianberg matrix based on the number of 
		// 	// // Arnoldi iterations
		// 	// for (int j = 0; j < restart_value + 1; j++) {
		// 	// 	if (j > i + 1) {
		// 	// 		free((*H)[j]);
		// 	// 	} else {
		// 	// 		(*H)[j] = realloc((*H)[j], sizeof(***H) * i);
		// 	// 	}
		// 	// }
		// 	// realloc((*H), sizeof(**H) * (i + 1));

		// 	// // reallocate residual based on the number of Arnoldi iterations
		// 	// g = realloc(g, sizeof(*g) * (i+1));

		// 	break;
		// }
	}

	// deallocate solver variables
	// mxFree(temp_basis_vector);
	// mxFree(c);
	// mxFree(s);
}

void copy_array(int n, double* array, double* copy) {
	for (int i = 0; i < n; i++) {
		copy[i] = array[i];
	}
}

void retrieve_array_from_matrix(int rows, int columns, double* matrix, 
	double* copy, int index, orientation_t orientation) {

	if (orientation == COLUMN) {
		int startingIndex = index * rows;
		for (int i = index * rows, j = 0; i < startingIndex + rows; i++, j++) {
			copy[j] = matrix[i];
		}

	} else if (orientation == ROW) {
		for (int i = index, j = 0; i < rows * columns; i = i + rows, j++) {
			copy[j] = matrix[i];
		}
	}
}
