#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <immintrin.h> // SSE4.1命令を使用するためのヘッダファイル
#include <pthread.h>  // マルチスレッドを使用するためのヘッダファイル
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>

// 行列の次元を設定
#define N 4
#define MATRIX_SIZE 4


// 行列掛け算関数
void matrix_multiply_direct(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double e[N][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            e[i][j]=0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                e[i][j] += A[i][k] * B[k][j];
            }
        }
    }

}

// 行列掛け算関数
void matrix_multiply(double A[MATRIX_SIZE][MATRIX_SIZE], double B[MATRIX_SIZE][MATRIX_SIZE], double C[MATRIX_SIZE][MATRIX_SIZE], int start_row, int end_row)
{
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            C[i][j] = 0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


// 行列の逆行列を計算する関数
void inverseMatrix(double A[N][N], double A_inv[N][N]) {
    int i, j, k;
    double temp;

    // 単位行列を初期化
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A_inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // ガウス・ジョルダン法による逆行列の計算
    for (k = 0; k < N; k++) {
        temp = A[k][k];
        for (j = 0; j < N; j++) {
            A[k][j] /= temp;
            A_inv[k][j] /= temp;
        }
        for (i = 0; i < N; i++) {
            if (i != k) {
                temp = A[i][k];
                for (j = 0; j < N; j++) {
                    A[i][j] -= A[k][j] * temp;
                    A_inv[i][j] -= A_inv[k][j] * temp;
                }
            }
        }
    }
}

int main() {
    double A[N][N];
    /*
     = {{1.0, 2.0, 3.0},
                      {4.0, 5.0, 6.0},
                      {7.0, 8.0, 10.0}};
    */
    double C[N][N];
    double A_inv[N][N];
    double AA[N][N];
    double e[N][N]={0};
    
    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++)
    A[i][j]=1.0+rand()%1129; //1.0+i*N+j;
    }
    A[2][2]=10.0;
    for(int i=0;i<N;i++){
    for(int j=0;j<N;j++)
    AA[i][j]=A[i][j];
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf,", A[i][j]);
        }
        printf("\n");
    }
    
    inverseMatrix(A, A_inv);


    printf("Inverse Matrix:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f\t", A_inv[i][j]);
        }
        printf("\n");
    }

    matrix_multiply_direct(AA,A_inv,e);
    printf("direct\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf,", e[i][j]);
        }
        printf("\n");
    }

printf("origin\n");
for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf,", AA[i][j]);
        }
        printf("\n");
    }
    

    printf("Inverse Matrix2:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f\t", A_inv[i][j]);
        }
        printf("\n");
    }
    //exit(1);

    int num_processes = 1; // 4つのプロセスを使用して行列掛け算を高速化

    // 各プロセスで一部の行を計算
    for (int i = 0; i < num_processes; i++)
    {
        pid_t pid = fork();

        if (pid == 0) // 子プロセスの処理
        {
            int start_row = i * (MATRIX_SIZE / num_processes);
            int end_row = (i + 1) * (MATRIX_SIZE / num_processes);

            matrix_multiply(AA, A_inv, e, start_row, end_row);

            // 結果を表示
            printf("Process %d: Rows %d to %d completed\n", i, start_row, end_row);

            exit(0); // 子プロセスを終了
        }
        else if (pid < 0)
        {
            // エラーハンドリング
            perror("fork");
            exit(1);
        }
    }

    //matrix_multiply_direct(AA,A_inv,e);
    // 親プロセスが子プロセスの終了を待つ
    for (int i = 0; i < num_processes; i++)
    {
        int status;
        wait(&status);
    }

    // 結果を表示（確認のため）
    printf("Result Matrix:\n");
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("%lf ", e[i][j]);
        }
        printf("\n");
    }

    return 0;
}

