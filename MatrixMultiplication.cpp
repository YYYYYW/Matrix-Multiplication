#include<stdio.h>
#include<stdlib.h>
#include<time.h>
static void mm_generate(float* matA,float* matB,float* matC,const int M,const int N,const int K)
{
	for (int i = 0; i < M;i++)
	{
		for (int j = 0; j < N;j++)
		{
			float sum = 0.0f;
			for (int k = 0; k < K;k++)
			{
				sum += matA[i*K + k] * matB[k*N + j];
			}
			matC[i*N + j] = sum;
		}
	}
}

static void showMatrix(float* C, int M, int N){
    printf("\n=============================\n");
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            printf("%f ",C[i*N + j]);
        }
        printf("\n");
    }
    printf("\n");
}

static void mm_strassen(float* matA, float* matB, float* matC, const int M, const int N, const int K)
{
	if ((M <= 2) || M%2 != 0 || N%2 != 0 || K%2!=0)
	{
		return mm_generate(matA, matB, matC, M, N, K);
	}

	int offset = 0;
	//M1 = (A11+A22)*(B11+B22)
	float* M1 = (float*) malloc((M/2) * (N/2) * sizeof(float));
	{
		//M1_0 = (A11+A22)
		float * M1_0 = (float*) malloc((M/2) * (K/2) * sizeof(float));
		offset = M*K / 2 + K / 2;
		for (int i = 0; i < M / 2; i++)
		{
			for (int j = 0; j < K/2; j++)
			{
				const int baseIdx = i*K + j;
				M1_0[i*K/2+j] = matA[baseIdx] + matA[baseIdx + offset];
			}
		}
		//M1_1 = (B11+B22)
		float* M1_1 = (float*) malloc((K/2) * (N/2) * sizeof(float));
		offset = K*N / 2 + N / 2;
		for (int i = 0; i < K / 2; i++)
		{
			for (int j = 0; j < N / 2; j++)
			{
				const int baseIdx = i*N + j;
				M1_1[i*N/2+j] = matB[baseIdx] + matB[baseIdx + offset];
			}
		}
		mm_strassen(&M1_0[0], &M1_1[0], &M1[0], M / 2, N / 2, K / 2);

		free(M1_0);         M1_0=NULL;
		free(M1_1);         M1_1=NULL;
	}

	//M2 = (A21+A22)*B11
	float* M2 = (float*) malloc((M/2) * (N/2) * sizeof(float));
	{
		//M2_0 = (A21+A22)
		float* M2_0 = (float*) malloc((M/2) * (K/2) * sizeof(float));
		offset = K / 2;
		for (int i = M / 2; i < M; i++)
		{
			for (int j = 0; j < K / 2; j++)
			{
				const int baseIdx = i*K + j;
				M2_0[(i-M/2)*K/2+j] = matA[baseIdx] + matA[baseIdx + offset];
			}
		}
		//M2_1 = B11
        float* M2_1 = (float*) malloc((K/2) * (N/2) * sizeof(float));
        for(int i = 0; i < K / 2; i++) {
            for(int j = 0; j < M / 2; j++){
                M2_1[i * N/2 + j] = matB[i * N + j];
            }
        }
		mm_strassen(&M2_0[0], &M2_1[0], &M2[0], M / 2, N / 2, K / 2);

		free(M2_0);         M2_0=NULL;
		free(M2_1);         M2_1=NULL;
	}

	//M3 = A11*(B12-B22)
	float* M3 = (float*) malloc((M/2) * (N/2) * sizeof(float));
	{
		//M3_0 = A11
		float* M3_0 = (float*) malloc((M/2) * (K/2) * sizeof(float));
		for(int i = 0; i < M / 2; i++){
            for(int j = 0; j < K / 2; j++){
                M3_0[i * K/2 + j] = matA[i * K + j];
            }
		}
		//M3_1 = (B12-B22)
		float* M3_1 = (float*) malloc((K/2) * (N/2) * sizeof(float));
		offset = K*N / 2;
		for (int i = 0; i < K/2; i++)
		{
			for (int j = N/2; j < N; j++)
			{
				const int baseIdx = i*N + j;
				M3_1[i*N/2+j-N/2] = matB[baseIdx] - matB[baseIdx + offset];
			}
		}
		mm_strassen(&M3_0[0], &M3_1[0], &M3[0], M / 2, N / 2, K / 2);

		free(M3_0);         M3_0=NULL;
		free(M3_1);         M3_1=NULL;
	}

	//M4 = A22*(B21-B11)
	float* M4 = (float*) malloc((M/2) * (N/2) * sizeof(float));
	{
		//M4_0 = A22
		float* M4_0 = (float*) malloc((M/2) * (K/2) * sizeof(float));
		for(int i = M / 2; i < M; i++){
            for(int j = K / 2; j < K; j++){
                M4_0[(i-M/2) * K/2 + j - K/2] = matA[i * K + j];
            }
		}
		//M4_1 = (B21-B11)
		float* M4_1 = (float*) malloc((K/2) * (N/2) * sizeof(float));
		offset = N*K/2;
		for (int i = 0; i < K / 2; i++)
		{
			for (int j = 0; j < N/2; j++)
			{
				const int baseIdx = i*N + j;
				M4_1[i*N/2 + j] = matB[baseIdx + offset] - matB[baseIdx];
			}
		}
		mm_strassen(&M4_0[0], &M4_1[0], &M4[0], M / 2, N / 2, K / 2);

		free(M4_0);         M4_0=NULL;
		free(M4_1);         M4_1=NULL;
	}

	//M5 = (A11+A12)*B22
	float* M5 = (float*) malloc((M/2) * (N/2) * sizeof(float));
	{
		//M5_0 = (A11+A12)
		float* M5_0 = (float*) malloc((M/2) * (K/2) * sizeof(float));
		offset = K / 2;
		for (int i = 0; i < M/2; i++)
		{
			for (int j = 0; j < K / 2; j++)
			{
				const int baseIdx = i*K + j;
				M5_0[i*K / 2 + j] = matA[baseIdx] + matA[baseIdx + offset];
			}
		}
		//M5_1 = B22
		float* M5_1 = (float*) malloc((K/2) * (N/2) * sizeof(float));
		offset = N*K/2 + N/2;
		for(int i = 0; i < K / 2; i++){
            for(int j = 0; j < N / 2; j++){
                M5_1[i * N/2 + j] = matB[i * N + j + offset];
            }
		}
		mm_strassen(&M5_0[0], &M5_1[0], &M5[0], M / 2, N / 2, K / 2);

		free(M5_0);         M5_0=NULL;
		free(M5_1);         M5_1=NULL;
	}

	//M6 = (A21-A11)*(B11+B12)
	float* M6 = (float*) malloc((M/2) * (N/2) * sizeof(float));
	{
		//M6_0 = (A21-A11)
		float* M6_0 = (float*) malloc((M/2) * (K/2) * sizeof(float));
		offset = K * M / 2;
		for (int i = 0; i < M / 2; i++)
		{
			for (int j = 0; j < K/2; j++)
			{
				const int baseIdx = i*K + j;
				M6_0[i*K/2+j] = matA[baseIdx + offset] - matA[baseIdx];
			}
		}
		//M6_1 = (B11+B12)
		float* M6_1 = (float*) malloc((K/2) * (N/2) * sizeof(float));
		offset = N / 2;
		for (int i = 0; i < K / 2; i++)
		{
			for (int j = 0; j < N/2; j++)
			{
				const int baseIdx = i*N + j;
				M6_1[i*N/2+j] = matB[baseIdx] + matB[baseIdx + offset];
			}
		}
		mm_strassen(&M6_0[0], &M6_1[0], &M6[0], M / 2, N / 2, K / 2);

		free(M6_0);         M6_0=NULL;
		free(M6_1);         M6_1=NULL;
	}

	//M7 = (A12-A22)*(B21+B22)
	float* M7 = (float*) malloc((M/2) * (N/2) * sizeof(float));
	{
		//M7_0 = (A12-A22)
		float* M7_0 = (float*) malloc((M/2) * (K/2) * sizeof(float));
		offset = M*K / 2;
		for (int i = 0; i < M / 2; i++)
		{
			for (int j = K/2; j < K; j++)
			{
				const int baseIdx = i*K + j;
				M7_0[i*K / 2 + j - K / 2] = matA[baseIdx] - matA[baseIdx + offset];
			}
		}
		//M7_1 = (B21+B22)
		float* M7_1 = (float*) malloc((K/2) * (N/2) * sizeof(float));
		offset = N / 2;
		for (int i = K/2; i < K; i++)
		{
			for (int j = 0; j < N / 2; j++)
			{
				const int baseIdx = i*N + j;
				M7_1[(i-K/2)*N / 2 + j] = matB[baseIdx] + matB[baseIdx + offset];
			}
		}
		mm_strassen(&M7_0[0], &M7_1[0], &M7[0], M / 2, N / 2, K / 2);

		free(M7_0);         M7_0=NULL;
		free(M7_1);         M7_1=NULL;
	}

	for (int i = 0; i < M / 2;i++)
	{
		for (int j = 0; j < N / 2;j++)
		{
			const int idx = i*N / 2 + j;
			//C11 = M1+M4-M5+M7
			matC[i*N + j] = M1[idx] + M4[idx] - M5[idx] + M7[idx];
			//C12 = M3+M5
			matC[i*N + j + N/2] = M3[idx] + M5[idx];
			//C21 = M2+M4
			matC[(i+M/2)*N + j] = M2[idx] + M4[idx];
			//C22 = M1-M2+M3+M6
			matC[(i+M/2)*N + j + N/2] = M1[idx] - M2[idx] + M3[idx] + M6[idx];
		}
	}
	free(M1);           M1=NULL;
	free(M2);           M2=NULL;
	free(M3);           M3=NULL;
	free(M4);           M4=NULL;
	free(M5);           M5=NULL;
	free(M6);           M6=NULL;
	free(M7);           M7=NULL;
}

static void mm_generate(float* matA, float* matB, float* matC, const int M, const int N, const int K,
                        const int strideA, const int strideB, const int strideC){
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            float sum = 0.0f;
            for(int k = 0; k < K; k++){
                sum += matA[i*strideA + k] * matB[k*strideB + j];
            }
            matC[i*strideC + j] = sum;
        }
    }
}

/*
 * matA M*K
 * matB K*N
 * matC M*N
 * matC = matA * matB
 * S1 = A21 + A22     T1 = B12 - B11
 * S2 = S1 - A11      T2 = B22 - T1
 * S3 = A11 - A21     T3 = B22 - B12
 * S4 = A12 - S2      T4 = T2 - B21
 * M1 = A11 * B11     U1 = M1 + M2
 * M2 = A12 * B21     U2 = M1 + M6
 * M3 = S4 * B22      U3 = U2 + M7
 * M4 = A22 * T4      U4 = U2 + M5
 * M5 = S1 * T1       U5 = U4 + M3
 * M6 = S2 * T2       U6 = U3 - U4
 * M7 = S3 * T3       U7 = U3 + M5
 * C11 = U1
 * C12 = U5
 * C21 = U6
 * C22 = U7
 */
static void mm_CoppersmithWinograd(float* matA, float* matB, float* matC, const int M, const int N, const int K,
                             const int strideA, const int strideB, const int strideC){
    if((M <= 2) || (M%2 != 0 || N%2 != 0 || K%2 != 0)){
        return mm_generate(matA, matB, matC, M, N, K, strideA, strideB, strideC);
    }

    float* S1 = (float*) malloc((M/2) * (K/2) * sizeof(float));
    float* S2 = (float*) malloc((M/2) * (K/2) * sizeof(float));
    float* S3 = (float*) malloc((M/2) * (K/2) * sizeof(float));
    float* S4 = (float*) malloc((M/2) * (K/2) * sizeof(float));
    {
        for(int i = 0; i < M/2; i++){
            for(int j = 0; j < K/2; j++){
                int idxA, offset, idxS = i * (K/2) + j;

                //S1     = A21 + A22
                idxA     = (i + (M/2)) * strideA + j;
                offset   = K/2;
                S1[idxS] = matA[idxA] + matA[idxA + offset];

                //S2     = S1 - A11
                idxA     = i * strideA + j;
                S2[idxS] = S1[idxS] - matA[idxA];

                //S3     = A11 - A21
                offset   = (M/2) * strideA;
                S3[idxS] = matA[idxA] - matA[idxA + offset];

                //S4     = A12 - S2
                idxA     = i * strideA + (K/2) + j;
                S4[idxS] = matA[idxA] - S2[idxS];
            }
        }
    }

    float* T1 = (float*) malloc((K/2) * (N/2) * sizeof(float));
    float* T2 = (float*) malloc((K/2) * (N/2) * sizeof(float));
    float* T3 = (float*) malloc((K/2) * (N/2) * sizeof(float));
    float* T4 = (float*) malloc((K/2) * (N/2) * sizeof(float));
    {
        for(int i = 0; i < K/2; i++){
            for(int j = 0; j < N/2; j++){
                int idxB, offset, idxT = i * (N/2) + j;

                //T1     = B12 - B11
                idxB     = i * strideB + j;
                offset   = (N/2);
                T1[idxT] = matB[idxB + offset] - matB[idxB];

                //T2     = B22 - T1
                idxB     = (i + (K/2)) * strideB + (N/2) + j;
                T2[idxT] = matB[idxB] - T1[idxT];

                //T3     = B22 - B12
                idxB     = i * strideB + (N/2) + j;
                offset   = ((K/2)) * strideB;
                T3[idxT] = matB[idxB + offset] - matB[idxB];

                //T4     = T2 - B21
                idxB     = (i + (K/2)) * strideB + j;
                T4[idxT] = T2[idxT] - matB[idxB];
            }
        }
    }

    //M1 = A11 * B11
    float* M1 = (float*) malloc((M/2) * (N/2) * sizeof(float));
    mm_CoppersmithWinograd(matA, matB, &M1[0], M/2, N/2, K/2, strideA, strideB, N/2);

    //M2 = A12 * B21
    float* M2 = (float*) malloc((M/2) * (N/2) * sizeof(float));
    mm_CoppersmithWinograd(&matA[K/2], &matB[(K/2)*strideB], &M2[0], M/2, N/2, K/2, strideA, strideB, N/2);

    //M3 = S4 * B22
    float* M3 = (float*) malloc((M/2) * (N/2) * sizeof(float));
    mm_CoppersmithWinograd(&S4[0], &matB[(K/2) * strideB + (N/2)], &M3[0], M/2, N/2, K/2, K/2, strideB, N/2);

    //M4 = A22 * T4
    float* M4 = (float*) malloc((M/2) * (N/2) * sizeof(float));
    mm_CoppersmithWinograd(&matA[(M/2) * strideA + (K/2)], &T4[0], &M4[0], M/2, N/2, K/2, strideA, N/2, N/2);

    //M5 = S1 * T1
    float* M5 = (float*) malloc((M/2) * (N/2) * sizeof(float));
    mm_CoppersmithWinograd(&S1[0], &T1[0], &M5[0], M/2, N/2, K/2, K/2, N/2, N/2);

    //M6 = S2 * T2
    float* M6 = (float*) malloc((M/2) * (N/2) * sizeof(float));
    mm_CoppersmithWinograd(&S2[0], &T2[0], &M6[0], M/2, N/2, K/2, K/2, N/2, N/2);

    //M7 = S3 * T3
    float* M7 = (float*) malloc((M/2) * (N/2) * sizeof(float));
    mm_CoppersmithWinograd(&S3[0], &T3[0], &M7[0], M/2, N/2, K/2, K/2, N/2, N/2);

    //C11 = U1 = M1 + M2
    //C12 = U5 = U4 + M3 = U2 + M5 + M3 = M1 + M6 + M5 + M3
    //C21 = U6 = U3 - M4 = U2 + M7 - M4 = M1 + M6 + M7 - M4
    //C22 = U7 = U3 + M5 = U2 + M7 + M5 = M1 + M6 + M7 + M5
    for(int i = 0; i < M/2; i++){
        for(int j = 0; j < N/2; j++){
            int idx = i * (N/2) + j;
            matC[i*strideC + j] = M1[idx] + M2[idx];
            matC[i*strideC + j + (N/2)] = M1[idx] + M6[idx] + M5[idx] + M3[idx];
            matC[(i+(M/2))*strideC + j] = M1[idx] + M6[idx] + M7[idx] - M4[idx];
            matC[(i+(M/2))*strideC + j + (N/2)] = M1[idx] + M6[idx] + M7[idx] + M5[idx];
        }
    }
    free(S1);           S1=NULL;
    free(S2);           S2=NULL;
    free(S3);           S3=NULL;
    free(S4);           S4=NULL;
    free(T1);           T1=NULL;
    free(T2);           T2=NULL;
    free(T3);           T3=NULL;
    free(T4);           T4=NULL;
    free(M1);           M1=NULL;
    free(M2);           M2=NULL;
    free(M3);           M3=NULL;
    free(M4);           M4=NULL;
    free(M5);           M5=NULL;
    free(M6);           M6=NULL;
    free(M7);           M7=NULL;
}

static void mm_CoppersmithWinograd(float* matA, float* matB, float* matC, const int M, const int N, const int K){
    mm_CoppersmithWinograd(matA, matB, matC, M, N, K, K, N, N);
}

void mm_test(int M, int N, int K, int rangeTop){
    unsigned seed = time(0);
    srand(seed);
    clock_t start,end;
    for(int i = 0; i < 10; i++){
        float * mA = (float*) malloc(M*K*sizeof(float));
        float * mB = (float*) malloc(K*N*sizeof(float));
        float * mC = (float*) malloc(M*N*sizeof(float));
        float * mD = (float*) malloc(M*N*sizeof(float));
        float * mE = (float*) malloc(M*N*sizeof(float));
        for(int j = 0; j < M*K; j++){
            mA[j] = rand() % rangeTop;
        }
        for(int j = 0; j < K*N; j++){
            mB[j] = rand() % rangeTop;
        }
        start = clock();
        mm_strassen(mA, mB, mC, M, N, K);
        end = clock();
        double endtime = (double) (end-start)/CLOCKS_PER_SEC;
        printf("Strassen%d time: %fms\n", i, endtime*1000);

        start = clock();
        mm_generate(mA, mB, mD, M, N, K);
        end = clock();
        endtime = (double) (end-start)/CLOCKS_PER_SEC;
        printf("Generate%d time: %fms\n", i, endtime*1000);

        start = clock();
        mm_CoppersmithWinograd(mA, mB, mE, M, N, K);
        end = clock();
        endtime = (double) (end-start)/CLOCKS_PER_SEC;
        printf("Winograd%d time: %fms\n", i, endtime*1000);

        for(int j = 0; j < M*N; j++){
            if(mC[j] != mD[j] || mC[j] != mD[j]){
                printf("========A========\n");
                showMatrix(mA, M, K);
                printf("========B========\n");
                showMatrix(mB, K, N);
                printf("========Strassen========\n");
                showMatrix(mC, M, N);
                printf("========Generate========\n");
                showMatrix(mD, M, N);
                printf("========Winograd========\n");
                showMatrix(mE, M, N);
                return ;
            }
        }
        printf("\n");
    }
}

int main(){;
    int M, N, K, rangeTop;
    M = 1000;
    N = 1000;
    K = 2000;
    rangeTop = 10;
    mm_test(M, N, K, rangeTop);
    return 0;
}
