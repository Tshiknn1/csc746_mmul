const char* dgemm_desc = "Basic implementation, three-loop dgemm.";

/*
 * This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in row-major format.
 * On exit, A and B maintain their input values.
 */
void square_dgemm(int n, double* A, double* B, double* C) 
{
   for (int i = 0; i < n; i++) { // rows
      for (int j = 0; j < n; j++) { // columns
         double square_sum = C[i*n + j];
         for (int k = 0; k < n; k++) {
            // C[i, j] = C[i, j] + A[i, n] * B[n, j];
            square_sum += A[i*n + k] * B[k*n + j];
         }
         C[i*n + j] = square_sum;
      }
   }
}
