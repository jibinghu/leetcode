import torch
import triton
import triton.language as tl

@triton.jit
def demo1_kernel(x_ptr, output_ptr):
    range = tl.arange(0, 8)
    mask = range < 5
    x = tl.load(x_ptr + range, mask, 0)
    # Store the results in output_ptr to retrieve later
    tl.store(output_ptr + range, x)

def run_demo1():
    print("Demo1 Output: ")
    input_tensor = torch.ones(8, dtype=torch.float32).cuda()  # Input tensor on GPU
    output_tensor = torch.zeros(8, dtype=torch.float32).cuda()  # Output tensor on GPU

    demo1_kernel[(1,)](input_tensor, output_tensor)  # Launch the kernel
    print("Kernel output:", output_tensor.cpu().numpy())  # Copy back to CPU for printing
    print("----------------------------------------------\n")

run_demo1()