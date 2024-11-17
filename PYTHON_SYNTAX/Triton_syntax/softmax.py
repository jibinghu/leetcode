import torch
import triton
import triton.language as tl


# Triton 内核：Softmax 前向计算
@triton.jit
def _softmax_fwd_kernel(
    output_ptr, output_stride,
    input_ptr, input_stride,
    n_cols,
    BLOCK_SIZE: tl.constexpr,
    num_warps: tl.constexpr
):
    # 获取线程块索引
    row_idx = tl.program_id(0)
    
    # 计算该线程块负责的行的起始地址
    input_row_start = input_ptr + row_idx * input_stride
    output_row_start = output_ptr + row_idx * output_stride

    # 定义本地块范围
    col_offsets = tl.arange(0, BLOCK_SIZE)
    
    # 仅处理 n_cols 范围内的列
    mask = col_offsets < n_cols

    # 加载输入行的内容到寄存器
    row = tl.load(input_row_start + col_offsets, mask=mask, other=-float('inf'))

    # 计算最大值（避免数值溢出）
    row_max = tl.max(row, axis=0)
    row = row - row_max

    # 指数计算
    numerator = tl.exp(row)
    denominator = tl.sum(numerator, axis=0)

    # 计算 Softmax 并写回输出
    result = numerator / denominator
    tl.store(output_row_start + col_offsets, result, mask=mask)


# 主函数：Softmax 的 Triton 实现
def softmax(x: torch.Tensor) -> torch.Tensor:
    """Triton 实现的 Softmax，只支持二维张量"""
    # 确保输入是二维张量
    rows, cols = x.shape
    assert x.dim() == 2, f"Expected 2D input, got {x.dim()}D input"

    # 计算 block_size，为大于等于 cols 的最小 2 的幂
    block_size = triton.next_power_of_2(cols)

    # 根据 block_size 动态调整线程数量 num_warps
    num_warps = 4  # 每个 warp 有 32 个线程
    if block_size > 2047:
        num_warps = 8
    if block_size > 4095:
        num_warps = 16

    # 定义网格大小，每个线程块（Block）处理一行数据
    grid = (rows,)

    # 创建与输入张量形状相同的空张量，用于存储输出
    sm_out = torch.empty_like(x)

    # 调用 Triton 内核
    _softmax_fwd_kernel[grid](
        sm_out, sm_out.stride(0),
        x, x.stride(0),
        cols,
        BLOCK_SIZE=block_size,
        num_warps=num_warps
    )

    return sm_out


# 测试 Softmax 函数
if __name__ == "__main__":
    # 创建一个随机输入张量
    x = torch.randn(8, 128, device='cuda')  # 8 行 128 列

    # 计算 Triton 加速的 Softmax
    sm_out = softmax(x)

    # 验证结果与 PyTorch 的 Softmax 实现是否一致
    # dim=1 表示沿着 第 1 维（也称为列维度）计算 Softmax，即对 每一行中的所有列 进行 Softmax 操作
    torch_out = torch.softmax(x, dim=1)

    # 对比误差
    print("Triton Softmax 与 PyTorch Softmax 的误差：", torch.allclose(sm_out, torch_out, atol=1e-6))