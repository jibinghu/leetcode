import triton
import triton.language as tl
import torch

@triton.jit
def demo3(z_ptr):
    range = tl.arange(0, 8)
    z = tl.store(z_ptr + range, 10, range < 5)


def run_demo3():
    print("Demo3 Output: ")
    z = torch.ones(4, 3).cuda()
    demo3[(1, 1, 1)](z)
    print(z.cpu().numpy())
          
run_demo3()