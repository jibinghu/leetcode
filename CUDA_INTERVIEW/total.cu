/***************************************************************************
 * 
 *  Copyright (C) 2024.10.8 ISCAS
 *  All rights reserved.
 * 
 *  File Name: Interview.cu
 *  Description: For ISCAS Use.
 * 
 *  This code, or any portion thereof, may not be reproduced, distributed,
 *  or transmitted in any form or by any means, including photocopying,
 *  recording, or other electronic or mechanical methods, without the prior
 *  written permission of the owner.
 * 
 *  Unauthorized use of this code or any part of it may result in legal
 *  action, and the owner reserves the right to pursue legal remedies to
 *  the fullest extent allowed by law.
 * 
 *  Contact Information: 544575367@qq.com
 * 
 ***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include "cuda.h"
#include "cuda_runtime.h"
#define threadsPerBlock 256
#include <iostream>
#include <vector>
using namespace std;
// ========================= reduce =====================
// naive reduce
__global__ void reduce (double *din, double *dout, int n) {
    int tid = threadIdx.x;
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    __shared__ double sdata [threadsPerBlock];
    if (id < n)
        sdata[tid] = din[id];
    __syncthreads();
    
    for (int s = 1; s < threadsPerBlock; s =  s << 1) {
        //s=1.就是0 2 4线程活着，把135写到自己
        if (tid % (2*s) == 0 && tid + s < threadsPerBlock) {
            sdata[tid] += sdata[tid+s];
        }
        __syncthreads();
    }
    //相比起来下面的更好，但这是逐渐变好的版本
    //一开始的时候 活着的线程数量更多，但后面活的少了，反而都是在不同warp里的
    //比如倒数第二次操作，0号去拿64，127去拿192的数，一个warp活跃的线程很少

    // for (int stride = threadsPerBlock / 2; stride >= 1; stride /= 2) {
    //     __syncthreads();
    //     if (tid < stride) {
    //         //前半线程把后半对应位置的数加到自己身上，
    //         //当stride小于32，也就是16开始，会bank冲突吗？
    //         //也没有，取时，只有前16线程取数，写时也是各自写
    //         sdata[tid] += sdata[tid+stride];
    //     }
    // }
    // __syncthreads();
    if (tid == 0) {
        dout[blockIdx.x] = sdata[0];
    }
    
}

// no warp divergence
__global__ void kernel1(double* arr, double* out, int N){
    __shared__ double s_data[threadsPerBlock];
    unsigned int tid = threadIdx.x;
    unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;  // tid号线程要负责的数组元素的位置
    if(i < N){
        s_data[tid] = arr[i];
    }
    __syncthreads();

    for(int s = 1; s < blockDim.x; s*=2){
        if(tid % (2*s) == 0 && i + s <N){     // 偶数线程work，
        //没看出好在哪？？？？？？？？？？？？？？？就是后面的虽然没用了但是还是工作
            s_data[tid] += s_data[tid + s];
        }
        __syncthreads();
    }

    if(tid == 0){
        out[blockIdx.x] = s_data[0];
    }
}
// dim3 gridSize ((N+255)/256)
// dim blockSize 256
// no bank conflict
//从这里开始看上面两个太low了
__global__ void kernel2(double* arr, double* out, int N){
    __shared__ double sdata[threadsPerBlock];
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int tid = threadIdx.x;
    sdata[tid] = arr[i];
    __syncthreads();
    // 优化1：这里后面32不用同步，可以使用寄存器通信
    // 优化2：可以让一个线程处理多个数据
    //for (int stride = threadsPerBlock / 2; stride >= 1; stride /= 2) {
    //这两句有细微的区别
    for (int s = threadsPerBlock / 2; s > 0; s = s >> 1) {

        if (tid < s) {
            sdata[tid] += sdata[tid+s];
        }
        __syncthreads();
    } 
    if (tid == 0) {
        out[blockIdx.x] = sdata[0];
    }
}


// 装逼版本1 
__global__ void reduce3(double * d_in, double * d_out, int total_num) {
    int tid = threadIdx.x;
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    __shared__ double sdata [threadsPerBlock];
    sdata[tid] = d_in[i];
    __syncthreads();
    for (int s = threadsPerBlock / 2; s > 32; s = s >> 1) { // 32 就结束
        if (tid < s && i + s < total_num) {
            sdata[tid] += sdata[tid+s];
        }
        __syncthreads();
    }

    if (tid < 32) { // 因为同一个warp所以不用__syncthreads同步
    //不做条件判断，每一步都全warp一起做 避免
        sdata[tid] += sdata[tid+32];
        sdata[tid] += sdata[tid+16];
        sdata[tid] += sdata[tid+8];
        sdata[tid] += sdata[tid+4];
        sdata[tid] += sdata[tid+2];
        sdata[tid] += sdata[tid+1];
        //每一步的取和写都没有bankconflict
    }
    if (tid == 0) {
        d_out[blockIdx.x] = sdata[0];
    }

}

////////////////////////------------=============================CURRENt
#define BLOCK_SIZE 256
// 装逼版本2 
template<typename T>
__global__ void reduce4(T * d_in, T * d_out, int N) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int tid = threadIdx.x;
    __shared__ T sdata [BLOCK_SIZE];
    sdata[tid] = d_in[i];
    __syncthreads();
    for (int s = BLOCK_SIZE / 2; s > 32; s = s >> 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid+s];
        }
        __syncthreads();
    }
    T sum = sdata[tid]; //得先取回到线程寄存器
    if (tid < 32) { // 使用shfl操作直接代替shared memory 
        //shfl只能是warp内使用，并且是寄存器通信
        //上一个版本仍是直接操作shared mem，那就不能用shfl
        
        sum += __shfl_down_sync(0xffffffff, sum, 16);
        sum += __shfl_down_sync(0xffffffff, sum, 8);
        sum += __shfl_down_sync(0xffffffff, sum, 4);
        sum += __shfl_down_sync(0xffffffff, sum, 2);
        sum += __shfl_down_sync(0xffffffff, sum, 1);
    }
    if (tid == 0) {
        d_out[blockDim.x] = sdata[0];
    }
}

// =================== reduce 2D ======================
// (M,N) -> (M,1)   首先 (M,N) -> (M,N/256) 再在CPU或者GPU上进行 (M,N/256) -> (M,1)
// gridSize&blockSize <<<(N/256,M), 256>>>
// 分block时，blockIdx.y指第几行，blockIdx.x指每行里，第几个256的块
__global__ void reduce2D_1(double * d_in, double * d_out, int N, int M) {
    int mid = blockIdx.y;
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int tid = threadIdx.x;
    __shared__ sdata [BLOCK_SIZE];
    sdata[tid] = d_in[mid][i];
    __syncthreads();
    for (int s = BLOCK_SIZE / 2; s >= 32; s = s >> 1) {
        if (tid < s && i + s < N) {
            sdata[tid] += sdata[tid+s];
        }
        __syncthreads();
    }
    double sum = sdata[tid];
    if (tid < 32) { // 使用shfl操作直接代替shared memory 
        //shfl只能是warp内使用因此上面没有办法用shfl
        sum += __shfl_down_sync(0xffffffff, sum, 16);
        sum += __shfl_down_sync(0xffffffff, sum, 8);
        sum += __shfl_down_sync(0xffffffff, sum, 4);
        sum += __shfl_down_sync(0xffffffff, sum, 2);
        sum += __shfl_down_sync(0xffffffff, sum, 1);
    }
    if (tid == 0) {
        d_out[mid][blockDim.x] = sdata[0];
    }
}

// (M,N) -> (1,N)   首先 (M,N) -> (M/256,N) 再在CPU或者GPU上进行 (M/256,N) -> (1,N)
// gridSize&blockSize <<<(M/256,N), 256>>>
__global__ void reduce2D_2(double * d_in, double * d_out, int N, int M) {
    int nid = blockIdx.y;
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int tid = threadIdx.x;
    __shared__ sdata [BLOCK_SIZE];
    sdata[tid] = d_in[i][nid];
    //不太好的是 直接就转换了下上一种方法的行列关系
    //同个warp取得数据是在列上连续
    //从global取数都没法连续，如果din是按列存的那还好
    __syncthreads();
    for (int s = BLOCK_SIZE / 2; s >= 32; s = s >> 1) {
        if (tid < s && i + s < M) {
            sdata[tid] += sdata[tid+s];
        }
        __syncthreads();
    }
    double sum = sdata[tid];
    if (tid < 32) { // 使用shfl操作直接代替shared memory 
        //shfl只能是warp内使用因此上面没有办法用shfl
        sum += __shfl_down_sync(0xffffffff, sum, 16);
        sum += __shfl_down_sync(0xffffffff, sum, 8);
        sum += __shfl_down_sync(0xffffffff, sum, 4);
        sum += __shfl_down_sync(0xffffffff, sum, 2);
        sum += __shfl_down_sync(0xffffffff, sum, 1);
    }
    if (tid == 0) {
        d_out[blockDim.x][nid] = sdata[0];
    }
}

// 当时来自阿里的追问
// Q1: 当 M 本身就小于 256 怎么办
// A1: 减少BLOCKSIZE,直到为32为止(BLOCKSIZE小于WARPSIZE明显不合适)
//想想 M行 所有的行reduce到一行，确实可能行很少
//我觉得上面每行连续的，都分到不同的block确实很不好

// Q2: 当M真的小于 32 怎么办
// A2: 修改代码 (目前只想到了暴力reduce)

// (M,N) -> (1,N) M < 32
// gridSize & blockSize <<<N/256, 256>>>
//为啥一定要256 我觉得32就更好，block起得多 多用sm硬件
//或者就每个warp负责一列，利用shfl可以比较好的把一列的数据进行reduce
__global__ void reduce2D_3(double * d_in, double * d_out, int N, int M) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    for (int j = 0; j < M; j++) {
        d_out[i] += d_in[j][i];
    }
}


// =================== findDiff ======================
// 在一群相同的数中找出唯一不同的数字，基本思路是前后比较，优化就是上shared memory，或者直接寄存器通信

//相同的数字异或会没掉，如果总共的数字数量是奇数，那相同数字有偶数个，一异或就变0，干脆就所有数字异或
//如果是偶数个 就在发现有俩数异或为0时记下来，最后thd0多异或他一次
//这里不用这么麻烦，或许下面这个方法真的还可以，但是操作有点蠢

__global__ void findDiff(int * d_in, int *d_out, int total_num) {
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    int front_id = (id + total_num - 1) % total_num;
    int next_id = (id + total_num + 1) % total_num;
    if (d_in[id] != d_in[front_id] && d_in[id] != d_in[next_id]) {
        d_out = d_in;
    }
}

// =================== DGEMM ======================
// naive GEMM
__global__ void DGEMM(double alpha, double beta, double * d_A, double * d_B, double * d_C, int m, int n, int k, int lda, int ldb, int ldc) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j = blockDim.y * blockIdx.y + threadIdx.y;
    double tmp = 0;
    //写太烂了 两个i是咋回事啊笑死
    for (int i = 0; i < k; i++) {
        tmp += d_A[i*lda+k] * d_B[k*ldb+j];
    }
    d_C[i*ldc+j] = alpha * tmp + beta * d_C[i*ldc+j];
}


// 使用shared memory且没有bank conflict的GEMM
// dim3 blcokSize (32,32)
// dim3 gridSize (M/32, N/32)
#define BLOCK_SIZE 32
__global__ void DGEMM3(double alpha, double beta, double * d_A, double * d_B, double * d_C, int m, int n, int k, int lda, int ldb, int ldc) {
    int bx = blockIdx.x;
    int by = blockIdx.y;
    int tx = threadIdx.x & 31;  //对32取余
    int ty = threadIdx.x >> 5;  //除32.同warp的 取得行也一样
    //block的每个线程都对shared mem上的tile进行1对1
    #define A(i,j) *((d_A) + (i) + (lda)*(j))
    #define B(i,j) *((d_B) + (i) + (ldb)*(j))
    #define C(i,j) *((d_C) + (i) + (ldc)*(j))
    #define sA(i,j) *((s_A) + (i) + (BLOCK_SIZE) *(j))  //j是行索引，i是列索引
    #define sB(i,j) *((s_B) + (i) + (BLOCK_SIZE) *(j))
    __shared__ double s_A[BLOCK_SIZE*BLOCK_SIZE];
    __shared__ double s_B[BLOCK_SIZE*BLOCK_SIZE];
    d_A += bx * BLOCK_SIZE; //在
    d_B += by * BLOCK_SIZE * ldb;
    d_C += bx * BLOCK_SIZE + by * BLOCK_SIZE * ldc;
    double tmp = 0;
    for (int inner_k = 0; inner_k < k; inner_k += BLOCK_SIZE) {
        sA(tx, ty) = A(tx, ty); 
        sB(ty, tx) = B(tx, ty); // 注意这里 sB(ty,tx) !!!!!!!!在global行存，到sharedmem column-major
        d_A += BLOCK_SIZE * lda;    //A矩阵在按行取？？？如果是列主序这里才是正确得
        d_B += BLOCK_SIZE;
        __syncthreads();
        for (int kk = 0; kk < BLOCK_SIZE; kk++) {
            tmp += sA(tx, kk) * sB(ty, kk); 
        }
        __syncthreads();//之后要取新得sA/B确保计算完成
    }
    C(tx, ty) = alpha * tmp + beta * C(tx, ty);
}


__global__ void DGEMM3(double alpha, double beta, double * d_A, double * d_B, double * d_C, int m, int n, int k, int lda, int ldb, int ldc) {
    int bx = blockIdx.x;
    int by = blockIdx.y;
    int tx = threadIdx.x & 31;  //对32取余
    int ty = threadIdx.x >> 5;  //除32.同warp的 取得行也一样
    //block的每个线程都对shared mem上的tile进行1对1
    #define A(i,j) *((d_A) + (i) + (lda)*(j))
    #define B(i,j) *((d_B) + (i) + (ldb)*(j))
    #define C(i,j) *((d_C) + (i) + (ldc)*(j))
    #define sA(i,j) *((s_A) + (i) + (BLOCK_SIZE) *(j))  //j是行索引，i是列索引
    #define sB(i,j) *((s_B) + (i) + (BLOCK_SIZE) *(j))
    __shared__ double s_A[BLOCK_SIZE*BLOCK_SIZE];
    __shared__ double s_B[BLOCK_SIZE*BLOCK_SIZE];
    d_A += bx * BLOCK_SIZE; //说明A的分块是沿着
    d_B += by * BLOCK_SIZE * ldb;
    d_C += bx * BLOCK_SIZE + by * BLOCK_SIZE * ldc;
    double tmp = 0;
    for (int inner_k = 0; inner_k < k; inner_k += BLOCK_SIZE) {
        sA(tx, ty) = A(tx, ty); 
        sB(ty, tx) = B(tx, ty); // 注意这里 sB(ty,tx) !!!!!!!!在global行存，到sharedmem column-major
        d_A += BLOCK_SIZE * lda;    //A矩阵在按行取？？？如果是列主序这里才是正确得
        d_B += BLOCK_SIZE;
        __syncthreads();
        for (int kk = 0; kk < BLOCK_SIZE; kk++) {
            tmp += sA(tx, kk) * sB(ty, kk); 
        }
        __syncthreads();//之后要取新得sA/B确保计算完成
    }
    C(tx, ty) = alpha * tmp + beta * C(tx, ty);
}

// =================== , ======================
// naive transpose
__global__ void transpose1(double * d_in, double * d_out, int N) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x < N && y < N) {
        d_out[x*N+y] = d_in[y*N+x]; 
    }
}
// 那就是合并访存，shared memory
// shared memory 合并访存
#define BLOCK_SIZE 32
__global__ void transpose2 (double * d_in , double * d_out, int M, int N) { ////！！！！
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int j = blockDim.y * blockIdx.y + threadIdx.y;
    __shared__ s_data[BLOCK_SIZE][BLOCK_SIZE+1];
    if (i < M && j < N ) {
        unsigned int index = i * N + j;
        s_data[threadIdx.y][threadIdx.x] = d_in[index]; //shared[内部j][内部i]=d_in[i][j]
    }
    __syncthreads();
    i = blockDim.y * blockIdx.y + threadIdx.x;  //现在各个block的索引都转置了，i还是用来定行号，但现在一行是M个数据了
    j = blockDim.x * blockIdx.x + threadIdx.y;  //现在用来指结果大矩阵里的列
    // 这个是连续的，32 32 连续，可以合并访存
    //写错了吧if (i < M && j < N ) {
    if (i < N && j < M ) {
        unsigned int index = i * M + j;
        d_out[index] = s_data[threadIdx.x][threadIdx.y];    //现在再按照shared mem里的各行的样子直接按位置搬到global
        //threadIdx在就是行！
    }
}

CURREN
// dim3 blockDim((M+31)/32, (N+31)/32);
// dim3 threadDim(32, 32)
// <<<blockDim, threadDim>>>
//blockIdx.x用来指引行
//从A里读时 是coalesced。然后按行放入sharedmem。取时按列取，再按行写回At，这样global的访问都是coalesced的
__global__ void transposeCoalesced(float *odata, const float *idata)
{
  __shared__ float tile[TILE_DIM][TILE_DIM];

  int x = blockIdx.x * TILE_DIM + threadIdx.x;  //指引列号
  int y = blockIdx.y * TILE_DIM + threadIdx.y;  //感觉y指引行号
  int width = gridDim.x * TILE_DIM;

  for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
     tile[threadIdx.y+j][threadIdx.x] = idata[(y+j)*width + x];

  __syncthreads();

  x = blockIdx.y * TILE_DIM + threadIdx.x;  //blockIdx.y在列方向上 定位块 threadIdx.x是在块内的列位置 transpose block offset
  y = blockIdx.x * TILE_DIM + threadIdx.y;  //转置后blockIdx.x * TILE_DIM用来指引行
    //为什么width不变啊啊啊啊
  for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
     odata[(y+j)*width + x] = tile[threadIdx.x][threadIdx.y + j];   //取得是shared mem里对应转置位置的数
     //唯一疑惑width应该不一样啊除非M==N
}

// ======================  GEMV =================================

// 这里M换成N，然后for循环去处理M。
// 同时如果M还是太大了w
template <unsigned int WarpSize>
__device__ __forceinline__ float warpReduceSum(float sum) {
    if (WarpSize >= 32)sum += __shfl_down_sync(0xffffffff, sum, 16); // 0-16, 1-17, 2-18, etc.
    if (WarpSize >= 16)sum += __shfl_down_sync(0xffffffff, sum, 8);// 0-8, 1-9, 2-10, etc.
    if (WarpSize >= 8)sum += __shfl_down_sync(0xffffffff, sum, 4);// 0-4, 1-5, 2-6, etc.
    if (WarpSize >= 4)sum += __shfl_down_sync(0xffffffff, sum, 2);// 0-2, 1-3, 4-6, 5-7, etc.
    if (WarpSize >= 2)sum += __shfl_down_sync(0xffffffff, sum, 1);// 0-1, 2-3, 4-5, etc.
    return sum;
}

// dim3 dimGrid(M/4);
// dim3 dimBlock(32,4); 
// 在问什么呢，什么变成 64 

// 这个哥们在干嘛
__global__ void Sgemv_v0( 
    float * __restrict__ A,
    float * __restrict__ x,
    float * __restrict__ y, 
    const int M,
    const int N) {
    // Block index
    int bx = blockIdx.x;    //在行方向分block，每个block一次负责blockDim.y行*blockDim.x列个的计算

    //按照grid和block的配置，一个block 4 warps，一个warp一行，就不需要shared mem了，warp内使用shfl

    // Thread index
    int tx = threadIdx.x;   // 不可能
    int ty = threadIdx.y;

    const int warp_size=32;  // warp_size 就是 32 !!!!!!!!!!!! 
    int laneId= tx;
    int current_row = blockDim.y * bx + ty;

    if(current_row < M){
        float ans=0;
        int kIteration = N/warp_size;   //n方向做多次

        for(int i=0; i< kIteration; i++){
            int current_col = i*warp_size + laneId;
            ans += A[current_row*N + current_col] * x[current_col];
        }
        //for(int current_col=laneId; current_col<N; current_col+warp_size)
        //我觉得换成这个也行
        ans += __shfl_down_sync(0xffffffff, ans, 16); // 0-16, 1-17, 2-18, etc.
        ans += __shfl_down_sync(0xffffffff, ans, 8);// 0-8, 1-9, 2-10, etc.
        ans += __shfl_down_sync(0xffffffff, ans, 4);// 0-4, 1-5, 2-6, etc.
        ans += __shfl_down_sync(0xffffffff, ans, 2);// 0-2, 1-3, 4-6, 5-7, etc.
        ans += __shfl_down_sync(0xffffffff, ans, 1);// 0-1, 2-3, 4-5, etc.
        if(laneId==0) y[current_row]=ans;
    }
}






// A: [32, 1024]

#define BLOCK_SIZE 256 /// 先写普通的一会再写shfl的
// =================== reduce 2D ======================
// (M,N) -> (M,1)   首先 (M,N) -> (M,N/256) 再在CPU或者GPU上进行 (M,N/256) -> (M,1)
// gridSize&blockSize <<<(N/256,M), 256>>>
// for 循环去做 (M,N/256) -> (M,1) 在GPU上，他是这个意思
// 
__global__ void reduce2D_1(double * d_in, double * d_out, int N, int M) {
    int mid = blockIdx.y;
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int tid = threadIdx.x;
    
    __shared__ sdata [BLOCK_SIZE];
    sdata[tid] = d_in[mid][i];
    __syncthreads();
    for (int s = BLOCK_SIZE / 2; s > 0; s = s >> 1) {
        if (tid < s && i + s < N) {
            sdata[tid] += sdata[tid+s];
        }
        __syncthreads();
    }
    if (tid == 0) {
        d_out[mid][blockDim.x] = sdata[0];
    }
    // double sum = sdata[tid];
    // if (tid < 32) { // 使用shfl操作直接代替shared memory 
    //     //shfl只能是warp内使用因此上面没有办法用shfl
    //     sum += __shfl_down_sync(0xffffffff, sum, 16);
    //     sum += __shfl_down_sync(0xffffffff, sum, 8);
    //     sum += __shfl_down_sync(0xffffffff, sum, 4);
    //     sum += __shfl_down_sync(0xffffffff, sum, 2);
    //     sum += __shfl_down_sync(0xffffffff, sum, 1);
    // }
    
}

// naive 版本， 
// 优化：1 用一个线程计算多个数据，
// 优化：2 rs开个数组，现在的reduce效率有点慢
__global__ void getPI(double *rs, int total_num) {
    int tid = threadIdx.x;
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int n = i * 2 + 1;
    // 把n用i替代
    double s = powf(-1, i);
    double tmp = s * 1.0 / n; //// n 
    // *rs += tmp; /// !!!!!!!!!!!!!!!!!!!!!!!!!
    atomicAdd(rs, tmp);
}


//看不懂问题是啥
__global__ void mergeArrays(const int* inputArray1, const int* inputArray2, int* outputArray, int size)
{
    // 这个是navie想法，问题就是atomicAdd性能太差
    
    // threadIdx 在这里，这个是输出的下标
    int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadIndex < size)
    {
        // index1 index2 初始化为0！！！！！！
        int value1 = inputArray1[index1];
        int value2 = inputArray2[index2];
        
        if (value1 < value2) {
            outputArray[threadIndex] = value1;
            atomicAdd(index1, 1);  // 使用atomicAdd
        } else {
            outputArray[threadIndex] = value2;
            atomicAdd(index2, 1);
        }
            
        
    }
}

// dim (32，,32)
// dim (M/32, N/32, ip2)
__global__ void matrixTranspose(float* input, float* output, int ip2, int n, int m) {
    // 定义共享内存
    __shared__ float sharedInput[BLOCK_SIZE][BLOCK_SIZE];

    // 计算当前线程的索引
    int bIndex = blockIdx.z;    //xi咋还有z？？？
    //还好，看起来是每一层的矩阵进行transpose，各层还是一样的
    int nIndex = blockIdx.y * blockDim.y + threadIdx.y;
    int mIndex = blockIdx.x * blockDim.x + threadIdx.x;

    // 计算输入矩阵和输出矩阵中的索引
    // 这里if去了，编译会报错！！！！！！！！！！！！！！！
    // 连续，
    int inputIndex = bIndex * n * m + nIndex * m + mIndex;  //nIndex是转之前的列号，看起来是按列存（每列m个连续存），mIndex行号
    int outputIndex = bIndex * n * m + mIndex * n + nIndex; //nIndex是转职后的行号，mIndex是列号，并且转置后，每列变成n个，这是符合的
    

    // 将输入矩阵元素加载到共享内存中
    //xi并且threadIdx.y原本指示列号，这里变shared里的行，是因为放入shared时，就是放转置后的样子，那就得按着sharedmem里的样子放进目标global
    if (nIndex < n && mIndex < m) {
        sharedInput[threadIdx.y][threadIdx.x] = input[inputIndex];
    }

    // 同步所有线程，确保共享内存加载完成
    __syncthreads();

    // 检查当前线程是否在矩阵的有效范围内
    if (nIndex < n && mIndex < m) {
        // 执行转置操作
        // no 不连续
        output[outputIndex] = sharedInput[threadIdx.x][threadIdx.y];    //xi我觉得不对，threadIdx.y在shreadmem里表示几列
        //写回output，threadIdx.y也应该是列号，
        //但是看了outputIndex= bIndex * n * m + mIndex * n + nIndex
        //nIndex = blockIdx.y * blockDim.y + threadIdx.y，说明threadIdx.y是global中的行号
    }
}


// 
__global__ void reduce2D_1(double * A, double * y, int lda, int n) {
    int mid = blockIdx.y;
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int tid = threadIdx.x;
    
    for (int t = 0; t < n / BLOCKSIZE; t++) {

    }
    __shared__ sdata [BLOCKSIZE];
    sdata[tid] = A[mid*lda+i];
    __syncthreads();
    for (int s = BLOCKSIZE / 2; s > 0; s = s >> 1) {
        if (tid < s && i + s < n) {
            sdata[tid] += sdata[tid+s];
        }
        __syncthreads();
    }
    // double sum = sdata[tid];
    // if (tid < 32) { // 使用shfl操作直接代替shared memory 
    //     //shfl只能是warp内使用因此上面没有办法用shfl
    //     sum += __shfl_down_sync(0xffffffff, sum, 16);
    //     sum += __shfl_down_sync(0xffffffff, sum, 8);
    //     sum += __shfl_down_sync(0xffffffff, sum, 4);
    //     sum += __shfl_down_sync(0xffffffff, sum, 2);
    //     sum += __shfl_down_sync(0xffffffff, sum, 1);
    // }
    if (tid == 0) {
        y[mid * lda + blockDim.x] = sdata[0];
    }
}

//  这个代码有bug，擦刚刚看出来
// 如果最大的数正好在最后一位，tid = 0 tid = 1 data 都会是那个最大值，先不要告诉她这个事情
// 没事，先这样说吧  这个也算比较聪明方法了
// 看他能不能看出来
// 第二小可能会有点问题
__global__ void findMinTwo(const int* arr, int* min1, int* min2) {
    int data = arr[threadIdx.x];
    int localMin1 = data, localMin2 = data;
    // 和那个没关系，就写我这个就可以，不是shared 的问题
    // 我们就寻找狭义的两个数
    // 这里后16线程是一直没有工作，他的意思可能是让
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {
        int other = __shfl_down_sync(0xFFFFFFFF, data, offset);
        if (other < data) {
            data = other;
        }
    }
    
//xixi有两个localMin来记录的话，不会让另个在同组里的数消失的吧

    4 5 1 3 5 6 7 8
    // 妈的不是说好了32个数组吗。
    // 
    if (threadIdx.x == 0) {
        *min1 = data;
    }

    data = arr[threadIdx.x];
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {
        int other = __shfl_down_sync(0xFFFFFFFF, data, offset);
        if (other < data && other != *min1) {
            data = other;
        }
    }

    if (threadIdx.x == 0) {
        *min2 = data;
    }

    


}
//当全为row major，黎课说mkn循环最佳（从数据局部性，A和B的数据近期访问的都是行数据
//MKN写法不仅可以优化访存，还可以优化程序的ILP。
//在MNK中，后一次K循环的乘积需要累加到前一次K循环的结果上，程序最内层循环存在迭代依赖。
//MKN的最内层循环不存在这样的迭代依赖，可以释放ILP。
void naive_row_major_sgemm_mkn(const float* A, const float* B, float* C, const int M,
    const int N, const int K) {
    int mi = 0;
    int ni = 0;
    int ki = 0;
    for (mi = 0; mi < M; mi ++) {
        for (ki = 0; ki < K; ki ++) {
            for (ni = 0; ni < N; ni ++) {
                C[mi * N + ni] += A[mi * K + ki] * B[ki * N + ni];
            }
        }
    }

}
// A [m,k]
// B [k,n]
// C [m,n]
void naive_row_major_sgemm(const float* A, const float* B, float* C, const int M,
    const int N, const int K) {
    for (int m = 0; m < M; ++m) {
        for (int n = 0; n < N; ++n) {
            T tmp = 0;
            for (int k = 0; k < K; ++k) {
                tmp += A[m * K + k] * B[k * N + n]; 
            }   
        } 
        C[m*N+n] = tmp;  //  // n在里面对于C矩阵是友好的
    }
}

// dim3 blockSize (32, 32)
// dim3 girdSize ((M+31)/32, (N+31)/32)
// 问问题是不是应该用 could 不用 can
// i can improve this implement  by storing A and B in shared memory
// because  both A and B will be reused 
// but for C, it can't be stored in shared memory
// 使用 reduce 进行优化，将每个A和B乘起来，然后reduce，他们乘积的和
__global__ void naive_row_major_sgemm(const float* A, const float* B, float* C, const int M,
    const int N, const int K) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int idx = threadIdx.x;
    int idy = threadIdx.y;
    __shared__ T sA [32*32];
    __shared__ T sB [32*32];
    // 其他都没有问题  是想重新写吗，这个是对的，我已经按照他的给你改好了
    if (i < M && j < N) {
        T tmp = 0;
        // 和右边是一样的！！
        // 淦，他看错了
        for (int innerk = 0; innerk < K; innerk += 32) { // 问一下idx是正确的吗，我感觉这里用Idx
                                                        // i 也是由idx组成的啊
            sA[idx*32 + idy] = A[i*K + inner_k + idy];  // 这里没有用合并访存，用jforM会更好
            sB[idx*32 + idy] = B[(innerk+idx) * N + j]; // megred memory access
            __syncthreads();
            for (int k = 0; k < 32; k++) {
                tmp += sA[idx*32+k] * B[k *32+ idy];  // 这里可以先将B进行transpose
            }
            __syncthreads();
        }
        C[i*N+j] = tmp;
    }
    
       
}

//xixi尝试优化
__global__ void better_naive_row_major_sgemm(const float* A, const float* B, float* C, const int M,
    const int N, const int K) {
    int i = blockIdx.y * blockDim.y;  //行块的起始位置
    int j = blockIdx.x * blockDim.x;  //列块的起始位置
    int idy = threadIdx.y;  //负责行
    int idx = threadIdx.x;  //负责列

    __shared__ T sA [32*32];
    __shared__ T sB [32*32];
    int tiles = (K+31)/32;
    // 其他都没有问题  是想重新写吗，这个是对的，我已经按照他的给你改好了
    if (i+idy < M && i+idx < N) {
        T tmp = 0;
        //k方向遍历获取宽度为32的块
        for (int k_tile = 0; k_tile < tiles; ++k_tile) { 
            if(k_tile*32+idx<N && k_tile*32+idy<N){
                sA[idy*32 + idx] = A[(i+idy)*K + k_tile*32+idx];  // xi，修改后合并访存了
                sB[idy*32 + idx] = B[(k_tile*32+idy)*N + i+ idx]; // 合并访存，而且也没有bank conflict
            }
            __syncthreads();
            for (int k = 0; k < 32; k++) {
                //一个warp的线程，每次会用相同的sA的数据，和sB同一行的不同的32个数，都没有bankconflict
                tmp += sA[idy*32+k] * sB[k *32+ idx];  
            }
            __syncthreads();
        }
        C[(i+idy)*N+j+idx] = tmp;
    }
    
}

/**
 * @brief 
 *  这个实验性的 用一个block计算 长度为 length数组的前缀和，Input为长度为n 的输入数组，只保证一个block计算前缀和
 *  output则是长度为n的Inpute 前缀和 pre_sum 前缀和 sum[i]= a[0]+a[1]+a[2]+..a[i]
 * @tparam index_t 
 * @tparam value_t 
 * @tparam warp_size  =32 
 * @param Input 
 * @param Output 
 * @param length 
 */
 template<typename index_t,typename value_t,int warp_size=32>
 __global__ void pre_sum_block(value_t * Input,value_t *Output,index_t length)
 {
     const int thid = blockDim.x*blockIdx.x + threadIdx.x; // 总的线程
     const int tx=threadIdx.x;
     const int wrapId = tx / warp_size;
     const int wraps =SDIV(blockDim.x,warp_size); // wraps<=32
     const int laneId = tx & (warp_size-1);// 取二进制最后五位，是 threadIdx对32取模的结果。
 
     if(thid>=length) return;
     // 越界
     value_t val = Input[thid]; // 每个线程的 负责一个数据，本地寄存器上
     __shared__ value_t pre_sum_block [32]; // 每个wrap的最后一个前缀和放在上面
     // const int iters = 
     // 计算 wrap内的前缀和
     #pragma unroll 5
     for(int delta=1;delta<warp_size;delta=delta*2) // 因为warp_size=32，否则应该是 delta< log2f(warp_size)
     {
          value_t temp=__shfl_up_sync(0xFFFFFFFF,val,delta,warp_size);
          if (laneId >=delta)
          //有分支
              val += temp;
         
     }
     // wrap是隐式同步的，限制每个wrap单独计算了前缀和
     if( laneId == warp_size-1)
     {
         // 一个wrap最后一个数
         pre_sum_block[wrapId]=val;
     }
     // 对shared memory的数求前缀和 ,wraps肯定是少于32的
     __syncthreads();// block内同步
 
 //给后面block的每个val加上前面的block的前缀和。又是一种前缀和
 //先用一个warp对pre_sum_block进行一个前缀和

     if(tx<warp_size) // 取第一个wrap对pre_sum_block计算
     {
         value_t warp_share_val = tx<wraps ?  pre_sum_block[tx] :0;
         #pragma unroll 5
         for(int delta=1;delta<warp_size;delta=delta*2) // 因为warp_size=32，否则应该是 delta< log2f(warp_size)
         {
             value_t temp=__shfl_up_sync(0xFFFFFFFF,warp_share_val,delta,warp_size);
             if (laneId >=delta)
                 warp_share_val += temp;
         }
 
         if(tx<wraps) 
             pre_sum_block[tx]= warp_share_val; // 每个wrap最后一个前缀和组成共享数组 的前缀和
 
     }
     __syncthreads();// block内同步，因为不同的warp要读share_memory
     if(wrapId>=1)  // 这里是 >=
     {
         //取wrap左边一个数
         val+=pre_sum_block[wrapId-1];
     }
     Output[thid]=val;
 }
 
 __global__ void scan(int *a,int *b,int equal_value, int N)
 {
    extern    __share__ int share_sum[];
    int tid=thread.x+blockIdx.x*blockDim.x;
    int temp1,temp=0;
    int i=0;
    int t_temp;
    int laneid=thread.x&0x1f,warpid=thread.x/warp_size;
    if((tid<N)&&(a[tid]==equal_value))
    {
      temp=1;
    }
    temp1=temp;//作为标记，用来标记是否写入 
  for(i=1;i<warp_size;i*=2)
  {
      t_temp = __shfl_up_sync(0xFFFFFFFF,temp,i,warp_size);
       if(laneid>=i)
       {
          temp+= t_temp;
       }               
   }
      //这里得出每个线程束的前缀和，且最后一个为最大
      if(laneid==(warp_size-1))
      {
          share_sum[warpid]=temp;
      }
      __sychthread();
      if(!tid)
      //每个block的0号线程把每个warp段的前缀和 加到后一个warp上
      {
           for(i=1;i<(N+blockDim.x-1)/warp_size;i++)
          {
              share_sum[i]=share_sum[i]+share_sum[i-1];
          }
      }
     __sychthread();
      if((laneid!=(warp_size-1))&&(warpid>0))
      { //为什么最后一个线程不要加，他也没有前面warp的前缀和啊
          temp+=share_sum[warpid-1];
      }
      __sychthread();
      if(temp1)
      {
          b[temp-1]=tid;
      }
}


//非CUDA=================================================================================================
//======================================================================================================


// ============================= CPU transpose ======================
struct Tensor {
	vector<int> data; // 张量具体数据
	vector<int> shapes; // 张量形状 {10, 10, 12, 30}  张量的初始perm {0,1,2,3} 
};
// 任意维度transpose函数   perm为转置后的排布 {2,1,0,3}
Tensor transpose(Tensor & d_in, const vector<int> & perm) {
	vector<int> t_data = d_in.data;
	vector<int> t_shapes = d_in.shapes;
	int dims = perm.size();
	Tensor d_out = d_in;
	int total_num = 1; // Tensor内数据元素的个数
    for (auto & i : t_shapes) {
		total_num *= i;
	}
    //统计元素总数

    //xi，先理解下tensor转置，比如原先10组12x30的数据，变成了12组10x30的数据，
    //变成将这10组,同一相对位置的30个数据拼在一起。[9,0,0]数据在新位置[0,9,0]
	
    for (int i = 0; i < total_num; i++) { // 依次遍历各个数据寻找其在d_out中的位置
		int next_i = 0;         // d_out中的位置
		int tmp_i = i;          // d_in中的位置
		vector<int> indexs(dims);   // 相对于shapes的indexs
		for (int t = dims-1; t >= 0; t--) { // 求indexs
			indexs[t] = tmp_i % t_shapes[t];    //对每一维取余
			tmp_i = tmp_i / t_shapes[t];    //除以 更小一维的dim，来看在更高一维里他的相对位置
		}
        //比如上述自己的例子t_shapes[10,12,30]，indexs就是记录获得他在各个维度上的位置

        int stride=1;
        for (int t = dims-1; t>= 0; t--) { // 确定目标的位置
			if (t != dims-1) {
				stride *= t_shapes[perm[t+1]];  //比如在计算新的第2维，那要先乘上新的第3维的stride（也就是dim)
			}
			next_i += indexs[perm[t]]*stride;
            //perm[3]=i,则在新维度的第三维的index 为旧维度的第i维的坐标，即indexs[i]
		}
        //xixi我改成上面！！！！
/*
		for (int t = dims-1; t>= 0; t--) { // 确定目标的位置
			if (t != dims-1) {
				next_i *= t_shapes[perm[t+1]];  //比如在计算新的第2维，那要先乘上新的第3维的stride（也就是dim)
			}
			next_i += indexs[perm[t]];

		}
*/
        //xi我感觉得是，新3维坐标+新2维坐标*新3维dim+新1维坐标*新2维dim*新3维dim
        //xi上面的循环感觉错了
		d_out.data[next_i] = d_in.data[i]; 
	}

    //xi完成struct中的成员shape填写
    for (int t = 0; t < dims; t++) {
        d_out.shapes[t] = t_shapes[perm[t]];
    }
	return d_out;	
}



//
#include 

int findMax(vector<int>& nums) {
    int len = nums.size();
    int left = 0;
    int right = nums.size() - 1;
    while (left < right) {
        int mid = left + (right - left) / 2;
        if (nums[mid] < nums[right]) {
            right = mid;
        }
        else {
            left = mid + 1;
        }
    }
    return (left-1)%len;
}


#include <iostream>
#include <vector>
using namespace std;
/// !!!!!!!!!!!!!!!!!!!!!! 函数写上面！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
int findMax(vector<int>& nums) {
    int len = nums.size();
    int left = 0;
    int right = nums.size() - 1;

    //xi看傻眼这啥？折半？？？
    while (left < right) {
        int mid = left + (right - left) / 2;
        if (nums[mid] < nums[right]) {
            right = mid; /////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        }
        else {
            left = mid + 1;   ////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        }
    }
    return (left-1+len)%len;
}
// idx 为5 他写错了
int main() {
	vector<int> vv = {4,5,6,7,8,9,1,2,3};
	cout << findMax(vv) << endl;
} 



// 

void delete_item (float * data, int len, int * len_ans) {
    //是说在len ans原地删除len个元素？
    if (len == 0) {
        return;
    }
    int fast = 1, slow = 1;
    while (fast < len) {
        if (data[fast] != data[fast - 1]) {
            //看样子是删除重复元素
            //那就双指针，把不重复的写在slow上
            data[slow] = data[fast];
            ++slow;
        }
        ++fast;
    }
    *len_ans = slow; /// *********************************************************
}


// n h w ip3 k_h k_w c_out d_w d_h s_h s_w bias   // 12个参数

// k_w_ k_h_ 是膨胀后的卷积核心
// k_w_ = d_w * (k_w-1) + 1
// k_h_ = d_h * (k_h-1) + 1


class Solution {
public:
    vector<int> inorderTraversal(TreeTreeNode* root) {
        vector<int> result;
        stack<TreeTreeNode*> st;
        if (root != NULL) st.push(root);
        while (!st.empty()) {
            TreeTreeNode* TreeNode = st.top();
            if (TreeNode != NULL) {
                st.pop(); // 将该节点弹出，避免重复操作，下面再将右中左节点添加到栈中
                if (TreeNode->right) st.push(TreeNode->right);  // 添加右节点（空节点不入栈）

                st.push(TreeNode);                          // 添加中节点
                st.push(NULL); // 中节点访问过，但是还没有处理，加入空节点做为标记。

                if (TreeNode->left) st.push(TreeNode->left);    // 添加左节点（空节点不入栈）
            } else { // 只有遇到空节点的时候，才将下一个节点放进结果集
                st.pop();           // 将空节点弹出
                TreeNode = st.top();    // 重新取出栈中元素
                st.pop();
                result.push_back(TreeNode->val); // 加入到结果集
            }
        }
        return result;
    }
};



// dp

#include <iostream>
#inlcude <vector>
using namespace std;
int lengthOfLIS(vector<int>& nums) {
    int maxRs = 0;
    int n = (int)nums.size();
    if (n == 0) {
        return 0;
    }
    vector<int> dp(n, 0);
    for (int i = 0; i < n; ++i) {
        dp[i] = 1;
        for (int j = 0; j < i; ++j) {
            if (nums[j] < nums[i]) {
                dp[i] = max(dp[i], dp[j] + 1);
            }
        }
        if (dp[i] > maxRs) {
            maxRs = dp[i];
        }
    }
    return maxRs;
}

int main()
{
	int n;
	int maxRs = 0;
	cin >> n;
	if (n == 0) return 0;
	vector<int> arr;
    arr.clear();///////
	for (int i = 0; i < n; i++) {
		int num;
		cin >> num;
		arr.emplace_back(num);
	}
	vector<int> dp(n, 0);
	for (int i = 0; i < n; ++i) {
        dp[i] = 1;
        for (int j = 0; j < i; ++j) {
            if (arr[j] < arr[i]) {
                dp[i] = max(dp[i], dp[j] + 1);
            }
        }
        if (dp[i] > maxRs) {
            maxRs = dp[i];
        }
    }
	cout << maxRs;
    return 0;
}



// sum 用 long  long ！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
pair<int,int> findOne(vector<int> & arr) {
    int n = arr.size() + 1;
    long long sum = (long long)n * (n+1)/2;
    long long mul = 1;
    for (int i = 1; i <= n; i++) {
        mul *= arr[i];
    }
    for (auto num : arr) {
        sum -= num;
        mul /= num;
    }
    
    // int ip1 = (sum-(int)(sqrt(double(sum*sum-4*mul))))/2;
	// int ip2 = (sum+(int)(sqrt(double(sum*sum-4*mul))))/2;
    int ip1 = (ans1-(int)(sqrt(double(ans1*ans1-4*ans2))))/2;
    int ip2 = (ans1+(int)(sqrt(double(ans1*ans1-4*ans2))))/2;
    return {ip1, ip2};
    
}

int main() {
    
    int n;
    cin >> n;
    for ()
}

class Solution {
public:
    vector<string> restoreIpAddresses(string s) {
        vector<string> ans;
        int n = s.length();
        //遍历IP的点可能的位置（第一个点）
        for(int i = 1; i < 4 && i < n - 2; i++){
            //第二个点的位置
            for(int j = i + 1; j < i + 4 && j < n - 1; j++){
                //第三个点的位置
                for(int k = j + 1; k < j + 4 && k < n; k++){
                    //最后一段剩余数字不能超过3
                    if(n - k >= 4)
                        continue;
                    // 从点的位置分段截取 分别为4个ip地址数据
                    //
                    string ip1 = s.substr(0, i);
                    string ip2 = s.substr(i, j - i);
                    string ip3 = s.substr(j, k - j);
                    string ip4 = s.substr(k);
                    //IP每个数字不大于255
                    if(stoi(ip1) > 255 || stoi(ip2) > 255 || stoi(ip3) > 255 || stoi(ip4) > 255)
                        continue;
                    //    排除前导0的情况
                    if((ip1.length() != 1 && ip1[0] == '0') || (ip2.length() != 1 && ip2[0] == '0') ||  (ip3.length() != 1 && ip3[0] == '0') || (ip4.length() != 1 && ip4[0] == '0'))
                        continue;
                    //组装IP地址
                    string temp = ip1 + "." + ip2 + "." + ip3 + "." + ip4;  /// !!!!! 双引号
                    ans.push_back(temp);
                }
            }
        }
        return ans; /// !!!!
    }
};


 template<typename T>
 class smart
 {
 private:
     T* _ptr;
     int* _count; //reference couting
 
 public:
     //构造函数
     smart(T* ptr = nullptr) :_ptr(ptr)
     {
         if (_ptr)
         {
             _count = new int(1);
         }
         else
         {
             _count = new int(0);
         }
     }
 
     //拷贝构造
     smart(const smart& ptr)
     {
         if (this != &ptr)
         {
             this->_ptr = ptr._ptr;
             this->_count = ptr._count;
 
             (*this->_count)++;
         }
     }
 
     //重载operator=
     smart& operator=(const smart & ptr)
     {
         if (this->_ptr == ptr._ptr)
         {
             return *this;
         }
         if (this->_ptr)
         {
             (*this->_count)--;
             if (*this->_count == 0)
             {
                 delete this->_ptr;
                 delete this->_count;
             }
         }
         this->_ptr = ptr._ptr;
         this->_count = ptr._count;
         (*this->_count)++;
         return *this;
     }
 
     //operator*重载
     T& operator*()
     {
         if (this->_ptr)
         {
             return *(this->_ptr);
         }
     }
 
     //operator->重载
     T* operator->()
     {
         if (this->_ptr)
         {
             return this->_ptr;
         }
     }
 
     //析构函数
     ~smart()
     {
         (*this->_count)--;
         if (*this->_count == 0)
         {
             delete this->_ptr;
             delete this->_count;
         }
     }
     //return reference couting
     int use_count()
     {
         return *this->_count;
     }
 };


#include<iostream>
#include<string.h>
using namespace std;

class String {
private:
	char* m_str;
public:
	// 无参构造
	String(const char* str = "") {
		// +1是为了包含\0
		int len = strlen(str) + 1;
		m_str = new char[len];
		strcpy_s(m_str, len, str);
	}
	// 拷贝构造
	String(const String& s) {
		int len = strlen(s.m_str) + 1;
		m_str = new char[len];
		strcpy_s(m_str, len, s.m_str);
	}
	// 析构
	~String() {
		if (m_str) {
			delete[] m_str;
			m_str = nullptr;
		}
	}
	// 赋值
	String& operator=(const String& s) {
		if (*m_str != *s.m_str) {
			if (m_str != nullptr) {
				delete[] m_str;
				m_str = nullptr;
			}
			int len = strlen(s.m_str) + 1;
			m_str = new char[len];
			strcpy_s(m_str, len, s.m_str);
		}
		return *this;
	}
};


#include<iostream>
#include<thread>
#include<mutex>
#include<condition_variable>
#include<queue>
using namespace std;


int main() {
	mutex mx;
	condition_variable cv;
	queue<int> q;
	const int capicity = 5;

	thread producer([&] {
		for (int i = 0; i < 10; i++) {
			unique_lock<mutex>lock(mx);
			cv.wait(lock, [&] {return q.size() <= capicity; });
			cout << "thread: " << this_thread::get_id() << "produce " << i << endl;
			q.push(i);
			cv.notify_all();
		}
		});

	thread consumer([&] {
		while (true) {
			unique_lock<mutex> lock(mx);
			cv.wait(lock, [&] {return !q.empty(); });
			cout << "thread: " << this_thread::get_id() << "consume " << q.front() << endl;
			q.pop();
			cv.notify_all();
		}
		});
	producer.join();
	consumer.join();
	return 0;
}



#include <iostream>
#include <queue>
// 这个好像他刚刚让你写一下
struct TreeNode {
    int value;
    TreeNode* left;
    TreeNode* right;
    TreeNode(int value):
        value(value), left(nullptr), right(nullptr) {}
};  // !!!!!!!!!!!!!!!  ;;;;;;;;;;;;;;;;;;;;
// 这个是可以的。直接写这个就ok
bool isCBT(TreeNode* head) {
    if (head == nullptr) {
        return true;
    }
    std::queue<TreeNode*> qtree; // 队列结构  
    qtree.push(head);
    TreeNode* tmp = nullptr;
    while (tmp = qtree.front()) {  // 将 tmp的左右结点依次入栈  做的是层级遍历！！！！！！！！ 逐层遍历的
        qtree.push(tmp->left);
        qtree.push(tmp->right);
        qtree.pop();   // 当前 tmp 结点出栈
    } 
    while(!qtree.empty()) {  // 上面是有null就停止了（检测到第一个null），但是剩下的null都还在队里
        if (qtree.front() != nullptr) {  // 按理来说，这边的应该都是null，如果有一个非null，说明这边的不是完全二叉树，因为按照这个遍历，只要出现第一个null，后面不可能有东西了
            return false;
        }
        qtree.pop();
    }
    return true;    // if pass the check, is CBT!
}

int main() {
    TreeNode* head1 = new TreeNode(1);
    head1->left = new TreeNode(2);
    head1->right = new TreeNode(3);
    head1->left->right = new TreeNode(4);
    head1->right->right = new TreeNode(5);

    std::cout << "==============CBT Test1==============\n";
    bool iscbt1 = isCBT(head1);
    std::cout << iscbt1 << std::endl;

    TreeNode* head2 = new TreeNode(1);
    head2->left = new TreeNode(2);
    head2->right = new TreeNode(3);
    head2->left->left = new TreeNode(4);
    head2->left->right = new TreeNode(5);
    head2->right->left = new TreeNode(6);

    std::cout << "==============CBT Test2==============\n";
    bool iscbt2 = isCBT(head2);
    std::cout << iscbt2 << std::endl;
    return 0;
}
// 网上没有原题
// 我目前想的是回溯
// 首先取第一位数字， 首先判断 是否在集合中有这个，如果有则先用（这里是一个回溯点，有可能这个数字是不对的） else 取比这个数更大的一个数字，这个情况下不要考虑回溯，其他数字都按照最小的来就可以
//
/// 从高到低，找到第一个数字在arr中不存在（而且还要比这个数字大），如果说在arr不存在，但是没有比这个数字大也是不可以的
int getNum(vector<int> arr, int m) {
    stack<int> ss;
    while (m > 0) {
        ss.push(m%10);
        m = m / 10;
    }
    while (!ss.empty()) { // 这里是stack，是从最高往最低的不是低到高
        int top = ss.top();
        int num = getMax(arr, top); // 获得大于等于top的最小数
        if (num == top) {
            // 这里应该写一个递归函数
            // 再去取下一个最大
        } else {
            // 这里直接就结束了剩下的都取最小的
        }
    }
}

// 

class Solution {
public:
    // 合并任意两个
    ListTreeNode* mergeTwoLists(ListTreeNode *h1, ListTreeNode *h2) {
        if ((!h1) || (!h2)) return h1 ? h1 : h2;
        ListTreeNode head, *tail = &head, *h1Ptr = h1, *h2Ptr = h2;
        while (h1Ptr && h2Ptr) {
            if (h1Ptr->val < h2Ptr->val) {
                tail->next = h1Ptr;
                h1Ptr = h1Ptr->next;
            } else {
                tail->next = h2Ptr;
                h2Ptr = h2Ptr->next;
            }
            tail = tail->next;
        }
        tail->next = (h1Ptr ? h1Ptr : h2Ptr);
        return head.next;
    }
    // 递归形式的reduce
    // zuobiyoubi 是这样的
    ListTreeNode* mergeReduce(vector <ListTreeNode*> &lists, int l, int r) {
        if (l == r) return lists[l];
        if (l > r) return nullptr;
        int mid = (l + r) >> 1;  /// l!!!!!!!!!!!!!!!!!!!!!!!
        return mergeTwoLists(mergeReduce(lists, l, mid), mergeReduce(lists, mid + 1, r));
    }
    //
    ListTreeNode* mergeKLists(vector<ListTreeNode*>& lists) {
        return mergeReduce(lists, 0, lists.size() - 1);
    }
};