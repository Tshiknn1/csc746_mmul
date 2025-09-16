#include <cstring>

const char* dgemm_desc = "Blocked dgemm.";

void square_dgemm_basic(int n, double* A, double* B, double* C) 
{
   for (int i = 0; i < n; i++) { // rows
      for (int j = 0; j < n; j++) { // columns
         for (int k = 0; k < n; k++) {
            // C[i, j] = C[i, j] + A[i, n] * B[n, j];
            C[i*n + j] += A[i*n + k] * B[k*n + j];
         }
      }
   }
}

void copy_block(double *dest, double *src, int n, int block_size, int row, int col) {
   for (int y = 0; y < block_size; y++) {
      std::memcpy(&dest[y * block_size],
         &src[(row * block_size + y) * n + col * block_size],
         block_size);
   }
}

void write_block(double *dest, double *src, int n, int block_size, int row, int col) {
   for (int y = 0; y < block_size, y++) {
      std::memcpy(&src[(row * block_size + y) * n + col * block_size],
         &dest[y * block_size],
         block_size);
   }
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in row-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm_blocked(int n, int block_size, double* A, double* B, double* C) 
{
   int Nb = n / block_size;
   int block_arr_size = block_size * block_size;
   double An[block_arr_size];
   double Bn[block_arr_size];
   double Cn[block_arr_size];
   for (int i = 0; i < Nb; i++) {
      for (int j = 0; j < Nb; j++) {
         copy_block(&Cn[0], C, n, block_size, i, j);

         for (int k = 0; k < Nb; k++) {
            copy_block(&An[0], A, n, block_size, i, k);
            copy_block(&Bn[0], B, n, block_size, k, j);
            square_dgemm_basic(block_size, &An[0], &Bn[0], &Cn[0]);
         }

         write_block(C, &Cn[0], n, block_size, i, j);
      }
   }
}
