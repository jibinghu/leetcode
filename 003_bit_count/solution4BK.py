class Solution:
    def countBits(self, n: int) -> list[int]:
        """
        最直观的做法是对从 000 到 nnn 的每个整数直接计算「一比特数」。每个 int\texttt{int}int 型的数都可以用 32 位二进制数表示，只要遍历其二进制表示的每一位即可得到 1 的数目。这里利用Brian Kernighan算法进行加速：对于任意整数 x，令 x=x & (x−1)x=x~\&~(x-1)x=x & (x−1)，该运算将 x 的二进制表示的最后一个 1 变成 0。因此，对 x 重复该操作，直到 x 变成 0，则操作次数即为 x 的「一比特数」。
        """
        def counts(x : int) -> int:
            nums = 0
            while(x > 0):
                x = x & (x - 1)
                nums+=1
            return nums
        return [counts(i) for i in range(0, n +1)]