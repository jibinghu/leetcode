cmake_minimum_required(VERSION 3.10)
project(ReduceV1Baseline LANGUAGES CXX CUDA)

# 设置 CUDA 标准
set(CMAKE_CUDA_STANDARD 11)

# 添加可执行文件
add_executable(reduce_v1_baseline reduce_v1_baseline.cu)

# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_60")