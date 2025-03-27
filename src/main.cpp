#include <iostream>
#include "slu_ddefs.h" // SuperLU header for double precision

int main() {
    // Define a 5x5 sparse matrix in CSC format (example)
    int n = 5;              // Matrix size
    int nnz = 12;           // Number of non-zeros
    double a[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0}; // Values
    int rowind[] = {0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4}; // Row indices
    int colptr[] = {0, 3, 6, 8, 10, 12}; // Column pointers

    // Right-hand side vector b
    double b[] = {17.0, 6.0, 13.0, 21.0, 21.0};

    // SuperLU variables
	SuperMatrix A, B, L, U;
    int *perm_r = (int*)malloc(n * sizeof(int));
    int *perm_c = (int*)malloc(n * sizeof(int));
    superlu_options_t options;
    SuperLUStat_t stat;
    int info;
	
    // Create matrix A in CSC format
    dCreate_CompCol_Matrix(&A, n, n, nnz, a, rowind, colptr, SLU_NC, SLU_D, SLU_GE);

    // Create dense matrix B for RHS
    dCreate_Dense_Matrix(&B, n, 1, b, n, SLU_DN, SLU_D, SLU_GE);

    // Set default options
    set_default_options(&options);
    options.ColPerm = NATURAL; // No column permutation for simplicity

    // Initialize statistics
    StatInit(&stat);

    // Factorize and solve
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

	if (info == 0) {
        double *x = ((double*)((DNformat*)B.Store)->nzval); // Corrected access
        std::cout << "Solution:\n";
        for (int i = 0; i < n; i++) {
            std::cout << x[i] << "\n";
        }
    } else {
        std::cout << "Error: info = " << info << "\n";
    }

    // Clean up
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);

    return 0;
}
