from typing import List
import sys

class Solution:
    def minSubArrayLen(self, s, nums:List[int]) -> int:
        sign = sys.maxsize
        for i in range(len(nums)):
            sums = 0
            for j in range(i, len(nums)):
                sums += nums[j]
                if sums >= s:
                    # current = j - i + 1
                    # if current < sign:
                    #     sign = current
                    sign = min(sign, j - i + 1)
                    break
        return sign if sign != sys.maxsize else 0
    
def main():
    solu = Solution()

    s = 7
    nums = [2, 3, 1, 2, 4, 3]
    print(solu.minSubArrayLen(s, nums))
    
if __name__ == "__main__":
    main()
