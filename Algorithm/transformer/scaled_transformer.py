import time
import torch
import torch.nn.functional as F

from sageattention import sageattn
F.scaled_dot_product_attention = sageattn

# 假设输入的 Query、Key 和 Value 张量维度为 (batch_size, num_heads, seq_len, embed_dim)
Q = torch.randn(4, 8, 128, 128).cuda()  # 例如: batch_size=4, num_heads=8, seq_len=64, embed_dim=64
K = torch.randn(4, 8, 128, 128).cuda()
V = torch.randn(4, 8, 128, 128).cuda()

start_time = time.time()
# 使用 scaled_dot_product_attention 计算 Attention 输出
attention_output = F.scaled_dot_product_attention(Q, K, V)

end_time = time.time()

print(f"消耗时间：{end_time - start_time}")
# print(attention_output.shape)  # 输出形状: (4, 8, 64, 64)
# print(attention_output)  # 输出形状: (4, 8, 64, 64)