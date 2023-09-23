#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>
#include <x86intrin.h> // SIMD命令を使用するためのヘッダファイル

#include "global-p.h"
#include "struct-p.h"

#define MATRIX_SIZE 8
#define SHM_KEY 1234


#define NUM_THREADS 4   // 使用するスレッドの数

// 行列を表示する関数
void print_matrix(int A[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%d ", A[i][j]);
        }
        printf("\n");
    }
}

// ガウスの消去法をSIMDで実行する関数
void gaussian_elimination_simd(double A[MATRIX_SIZE][MATRIX_SIZE]) {
    for (int k = 0; k < MATRIX_SIZE; k++) {
        for (int i = k + 1; i < MATRIX_SIZE; i++) {
            __m128d factor = _mm_set1_pd(A[i][k] / A[k][k]);
            for (int j = k; j < MATRIX_SIZE; j += 2) {
                __m128d row1 = _mm_loadu_pd(&A[i][j]);
                __m128d row2 = _mm_loadu_pd(&A[k][j]);
                row1 = _mm_sub_pd(row1, _mm_mul_pd(factor, row2));
                _mm_storeu_pd(&A[i][j], row1);
            }
        }
    }
}


// スレッドごとに処理する行の範囲を指定する構造体
typedef struct {
    int (*matrix)[MATRIX_SIZE];
    int start_row;
    int end_row;
} ThreadData;

// ガウスの消去法を実行する関数（マルチスレッド用）
void *gaussian_elimination_thread(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    int (*A)[MATRIX_SIZE] = data->matrix;
    int start_row = data->start_row;
    int end_row = data->end_row;

    // ガウスの消去法をSIMDで実行
    gaussian_elimination_simd(A);

    pthread_exit(NULL);
}

int maina_gausu(MTX A){
    //int A[MATRIX_SIZE][MATRIX_SIZE]; // 入力行列
    // 行列 A を初期化（例：ランダムな値を設定）
    /*
	for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A.x[i][j] = rand() % 100; // 0から99までのランダムな整数
        }
    }
	*/
    printf("Input Matrix:\n");
    print_matrix(A.x);

    // マルチスレッドのセットアップ
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];
    int rows_per_thread = MATRIX_SIZE / NUM_THREADS;

    // 各スレッドに行列 A の一部を処理させる
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].matrix = A.x;
        thread_data[i].start_row = i * rows_per_thread;
        thread_data[i].end_row = (i + 1) * rows_per_thread;
        pthread_create(&threads[i], NULL, gaussian_elimination_thread, &thread_data[i]);
    }

    // 各スレッドの終了を待つ
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("Result Matrix:\n");
    print_matrix(A.x);
}

// 行列の掛け算関数
void matrix_multiply(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double *C, int start_row, int end_row) {
    for (int i = start_row; i < end_row; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            double sum = 0.0;
            for (int k = 0; k < MATRIX_SIZE; k++) {
                sum += A[i][k] * B[k][j];
            }
            C[i*MATRIX_SIZE+j] = sum;
        }
    }
}


// 行列の逆行列を計算する関数 (SSE4.1 SIMD版)
void matrix_inverse_simd(int A[MATRIX_SIZE][MATRIX_SIZE], double result[MATRIX_SIZE][MATRIX_SIZE]) {
    int pivot[MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; i++) {
        pivot[i] = -1;
    }

    __m128d one = _mm_set1_pd(1.0);

    for (int col = 0; col < MATRIX_SIZE; col++) {
        int pivot_row = -1;
        double max_value = 0.0;

        for (int row = 0; row < MATRIX_SIZE; row++) {
            if (pivot[row] != -1) continue;

            double val = fabs(A[row][col]);
            if (val > max_value) {
                max_value = val;
                pivot_row = row;
            }
        }

        if (pivot_row == -1) {
            fprintf(stderr, "Matrix is singular.\n");
            return;
        }

        pivot[pivot_row] = col;

        // Scale the pivot row
        double pivot_value = A[pivot_row][col];
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[pivot_row][j] /= pivot_value;
            result[pivot_row][j] = A[pivot_row][j];
        }

        // Eliminate non-zero entries below the pivot
        for (int row = 0; row < MATRIX_SIZE; row++) {
            if (row == pivot_row) continue;

            __m128d scale = _mm_set1_pd(A[row][col]);
            for (int j = 0; j < MATRIX_SIZE; j += 2) {
                __m128d row_pivot = _mm_loadu_pd(result[pivot_row] + j);
                __m128d scaled = _mm_mul_pd(scale, row_pivot);
                __m128d row_target = _mm_loadu_pd(result[row] + j);
                row_target = _mm_sub_pd(row_target, scaled);
                _mm_storeu_pd(result[row] + j, row_target);
            }
        }
    }
}

// 行列の逆行列を計算する関数
void inverseMatrix(double A[MATRIX_SIZE][MATRIX_SIZE], double A_inv[MATRIX_SIZE][MATRIX_SIZE]) {
    int i, j, k;
    double temp;

    // 単位行列を初期化
    for (i = 0; i < MATRIX_SIZE; i++) {
        for (j = 0; j < MATRIX_SIZE; j++) {
            A_inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // ガウス・ジョルダン法による逆行列の計算
    for (k = 0; k < MATRIX_SIZE; k++) {
        temp = A[k][k];
        for (j = 0; j < MATRIX_SIZE; j++) {
            A[k][j] /= temp;
            A_inv[k][j] /= temp;
        }
        for (i = 0; i < MATRIX_SIZE; i++) {
            if (i != k) {
                temp = A[i][k];
                for (j = 0; j < MATRIX_SIZE; j++) {
                    A[i][j] -= A[k][j] * temp;
                    A_inv[i][j] -= A_inv[k][j] * temp;
                }
            }
        }
    }
}

// 行列掛け算関数（SIMD版、SSE4.1）
void matrix_multiply_simd(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], int *C, int start_row, int end_row)
{
    
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j += 4) // 4つの要素を同時に処理するために4で割る
        {
            __m128d sum = _mm_setzero_pd(); // ゼロで初期化された128ビット整数ベクトルを作成

            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                __m128d a = _mm_set1_pd(A[i][k]); // 行列Aの要素をベクトル化
                __m128d b = _mm_loadu_pd((__m128d *)&B[k][j]); // 行列Bの要素をベクトル化

                sum = _mm_add_pd(sum, _mm_mul_pd(a, b)); // ベクトル同士の掛け算と足し算
            }

            // ベクトルの要素を合算して結果を格納
            double result[4];
            _mm_storeu_pd((__m128d *)result, sum);
            for (int k = 0; k < 4; k++)
            {
                C[i * MATRIX_SIZE + j + k] = result[k];
				C[i*MATRIX_SIZE+j+k]%=N;
            }
        }
    }
}



int matmul_simd(double matrixA[MATRIX_SIZE][MATRIX_SIZE],double matrixB[MATRIX_SIZE][MATRIX_SIZE]) {
    double resultMatrix[MATRIX_SIZE][MATRIX_SIZE] = {0};

    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            __m128d sum = _mm_setzero_pd();

            for (int k = 0; k < MATRIX_SIZE; k++) {
                __m128d a = _mm_set1_pd(matrixA[i][k]);
                __m128d b = _mm_loadu_pd(&matrixB[k][j]);
                sum = _mm_add_pd(sum, _mm_mul_pd(a, b));
            }

            _mm_storeu_pd(&resultMatrix[i][j], sum);
        }
    }

    // 結果を表示
    printf("Result Matrix:\n");
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%lf ", resultMatrix[i][j]);
        }
        printf("\n");
    }

    return 0;
}

//#include <immintrin.h>
//#include <stdio.h>

//#define MATRIX_SIZE 2

// 行列の逆行列を計算する関数
void inverseMatrix_simd(double A[MATRIX_SIZE][MATRIX_SIZE], double A_inv[MATRIX_SIZE][MATRIX_SIZE]) {
    int i, j, k;
    double temp;

    // 単位行列を初期化
    for (i = 0; i < MATRIX_SIZE; i++) {
        for (j = 0; j < MATRIX_SIZE; j++) {
            A_inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // ガウス・ジョルダン法による逆行列の計算
    for (k = 0; k < MATRIX_SIZE; k++) {
        temp = A[k][k];
        for (j = 0; j < MATRIX_SIZE; j++) {
            A[k][j] /= temp;
            A_inv[k][j] /= temp;
        }
        for (i = 0; i < MATRIX_SIZE; i++) {
            if (i != k) {
                temp = A[i][k];
                for (j = 0; j < MATRIX_SIZE; j++) {
                    A[i][j] -= A[k][j] * temp;
                    A_inv[i][j] -= A_inv[k][j] * temp;
                }
            }
        }
    }
}



int main() {
    double A[MATRIX_SIZE][MATRIX_SIZE];
    double A_inv[MATRIX_SIZE][MATRIX_SIZE];
    double C[MATRIX_SIZE][MATRIX_SIZE];
    double AA[MATRIX_SIZE][MATRIX_SIZE];

    // 行列 A を初期化
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i][j] = 1.0 + rand() % 7;
            printf("%lf,",A[i][j]);
            AA[i][j]=A[i][j];
        }
        printf("\n");
    }

    // 行列 A のコピーを作成
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            AA[i][j] = A[i][j];
        }
    }

    inverseMatrix_simd(A,A_inv);
    for(int i=0;i<MATRIX_SIZE;i++){
    printf(":%d ",i);
    for(int j=0;j<MATRIX_SIZE;j++)
    printf("%lf,",A_inv[i][j]);
    printf("\n");
    }
    matmul_simd(AA,A_inv);
    exit(1);

    // マルチプロセスで行列掛け算を並列化
    int num_processes = 8;
    int rows_per_process = MATRIX_SIZE / num_processes;

    int shmid = shmget(SHM_KEY, sizeof(double) * MATRIX_SIZE * MATRIX_SIZE, IPC_CREAT | 0666);
    if (shmid == -1) {
        perror("shmget");
        exit(1);
    }

    double *shared_C = (double *)shmat(shmid, NULL, 0);
    if (shared_C == (double *)-1) {
        perror("shmat");
        exit(1);
    }

    // 各プロセスで一部の行を計算
    for (int i = 0; i < num_processes; i++) {
        pid_t pid = fork();

        if (pid == 0) {
            int start_row = i * rows_per_process;
            int end_row = (i + 1) * rows_per_process;


            matrix_multiply_simd(AA,A_inv,shared_C,start_row,end_row);
            //matrix_multiply(AA, A_inv, shared_C, start_row, end_row);

            // 結果を表示
            printf("Process %d: Rows %d to %d completed\n", i, start_row, end_row);

            exit(0);
        } else if (pid < 0) {
            perror("fork");
            exit(1);
        }
    }

    // 親プロセスが子プロセスの終了を待つ
    for (int i = 0; i < num_processes; i++) {
        int status;
        wait(&status);
    }

    double D[MATRIX_SIZE][MATRIX_SIZE];
    // 結果を表示
    printf("Result Matrix:\n");
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%lf ", shared_C[i * MATRIX_SIZE + j]);
        }
        printf("\n");
    }

    // 共有メモリを解放
    if (shmdt(shared_C) == -1) {
        perror("shmdt");
        exit(1);
    }
    if (shmctl(shmid, IPC_RMID, NULL) == -1) {
        perror("shmctl");
        exit(1);
    }

    return 0;
}
