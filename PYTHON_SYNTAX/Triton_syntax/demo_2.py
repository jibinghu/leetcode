import torch
import triton
import triton.language as tl

@triton.jit
def demo2(x_ptr, x_out_ptr):
    i_range = tl.arange(0, 8)[:, None]  # shape (8, 1)
    j_range = tl.arange(0, 4)[None, :]  # shape (1, 4)
    range = i_range * 4 + j_range  # shape (8, 4)

    # 读取输入数据，条件是掩码范围
    mask = (i_range < 4) & (j_range < 3)  # 掩码：有效范围在前 4 行和前 3 列
    x = tl.load(x_ptr + range, mask, 0)   # 从输入指针读取数据，超出范围填充 0

    # 写入输出数据
    tl.store(x_out_ptr + range, x)  # 将读取到的数据写入到输出张量

def run_demo2():
    print("Demo2 Output: ")
    input_tensor = torch.ones(4, 4).cuda()  # 输入张量 (4, 4)，全为 1
    out = torch.empty(8, 4).cuda()          # 输出张量 (8, 4)，初始化为空
    demo2[(1, 1, 1)](input_tensor, out)    # 启动 Triton 内核
    print(out.cpu().numpy())               # 将结果从 GPU 移动到 CPU 并打印

run_demo2()