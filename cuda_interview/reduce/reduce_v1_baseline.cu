/*
* 总述：
* 通过不断地将多个元素合并成一个元素，从而减少数据量。场景：计算总和、最小值、最大值等。
* 在 GPU 中，通常利用树状结构进行 reduce 操作； 
* 假设给定一个长度为 N 的数组，需要计算所有元素之和。首先将数组分为 m 小份，
* 第一阶段中开启 m 个 block 计算出 m 个 reduce 值；然后在第二阶段用一个 block 将 m 个值再次 reduce
* Adapted by https://github.com/Liu-xiandong/How_to_optimize_in_GPU
*/

#include <cuda_runtime.h>
#include <cuda.h>
#include <stdio.h>
#include <math.h>

#define THREAD_PER_BLOCK 1024 // A100 系列的 THREADS 数量限制；

// 最基础版本，仅用于个人学习；
__global__ void reduce (double *d_in, double *d_out, int n) {
    int tid = threadIdx.x;
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    __shared__ double sdata [THREAD_PER_BLOCK];

    if (id < n)
        sdata[tid] = d_in[id];
    __syncthreads();
    // s 作为递归扩张系数
    for (int s=1; s<THREAD_PER_BLOCK; s=s<<1) {
        // 这里的 tid+s<THREAD_PER_BLOCK 事实上有些鸡肋：
        // 考虑临界值的情况下， s 再扩大一倍的话必然不满足 s<THREAD_PER_BLOCK 条件；
        if (!(tid%(2*s)) && tid+s<THREAD_PER_BLOCK) {
            sdata[tid] += sdata[tid+s];
        }
        __syncthreads();
    }
    // 注意这里将每个 BLOCK 中归约后的值保存到 d_out 中，也就是 d_out 的大小一定是 N/THREAD_PER_BLOCK
    if (tid == 0) {
        d_out[blockIdx.x] = sdata[0];
    }
}

// 定义 check 函数
bool check(double *out, double *res, int N){
    for(int i=0; i<N/THREAD_PER_BLOCK; i++){
        if(fabs(out[i]-res[i])>1e-6)
            return false;
    }
    return true;
}

// 定义 main 函数
int main(){
    const int N = 1024;
    // 为 input 分配空间
    double *in = (double *)malloc(N*sizeof(double));
    double *d_in;
    cudaMalloc((void **)&d_in,N*sizeof(double)) ;
    // 为 output 数组分配空间
    double *out = (double *)malloc(N/THREAD_PER_BLOCK*sizeof(double));
    double *d_out;
    cudaMalloc((void **)&d_out, N/THREAD_PER_BLOCK*sizeof(double));
    // 初始化原数组
    for(int i=0; i<N; i++){
        in[i] = 1.0;
    }
    // 将原数组复制到设备
    cudaMemcpy(d_in,in,N*sizeof(double),cudaMemcpyHostToDevice);
    // 创建 grid/block，实际运行
    dim3 grid(N/THREAD_PER_BLOCK, 1);
    dim3 block(THREAD_PER_BLOCK, 1);
    reduce<<<grid,block>>>(d_in, d_out, N);
    // 复制 output 到主机
    cudaMemcpy(out,d_out,N/THREAD_PER_BLOCK*sizeof(double),cudaMemcpyDeviceToHost);
    // 创建验证矩阵并验证
    double *res = (double *)malloc(N/THREAD_PER_BLOCK*sizeof(double));
    for(int i=0;i<N/THREAD_PER_BLOCK;i++){
        double cur = 0.0;
        for(int j=0;j<THREAD_PER_BLOCK;j++){
            cur += in[i*THREAD_PER_BLOCK+j];
        }
        res[i] = cur;
    }

    if(check(out,res,N))printf("the ans is right\n");
    else{
        printf("the ans is wrong\n");
        for(int i=0;i<N/THREAD_PER_BLOCK;i++){
            printf("%lf ",out[i]);
        }
        printf("\n");
    }

    free(in);
    free(out);
    free(res);
    cudaFree(d_in);
    cudaFree(d_out);
}