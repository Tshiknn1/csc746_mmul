#include <cstring>
#include <iostream>

const char* dgemm_desc = "Blocked dgemm.";

void square_dgemm_basic(int n, double* A, double* B, double* C) 
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

void copy_block(double *dest, double *src, int n, int block_size) {
   for (int y = 0; y < block_size; y++) {
      std::memcpy(&dest[y * block_size],
         &src[y * n],
         block_size);
   }
}

void write_block(double *dest, double *src, int n, int block_size) {
   for (int y = 0; y < block_size; y++) {
      std::memcpy(&src[y * block_size],
         &dest[y * n],
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

   double* An = new double[block_arr_size];
   double* Bn = new double[block_arr_size];
   double* Cn = new double[block_arr_size];

   for (int i = 0; i < Nb; i++) {   // row
      for (int j = 0; j < Nb; j++) {   // col
         int Cpos = i * n * block_size + j * block_size;
         copy_block(&Cn[0], &C[Cpos], n, block_size);

         for (int k = 0; k < Nb; k++) {
            copy_block(An, &A[k * n * block_size + j * block_size], n, block_size);
            copy_block(Bn, &B[i * n * block_size + k * block_size], n, block_size);

            square_dgemm_basic(block_size, An, Bn, Cn);
         }

         write_block(&C[Cpos], &Cn[0], n, block_size);
      }
   }

   delete[] An;
   delete[] Bn;
   delete[] Cn;
}
