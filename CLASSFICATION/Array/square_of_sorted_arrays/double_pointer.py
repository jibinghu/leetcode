from typing import List
import unittest

class Solution:
    def sortedSquares(self, nums:List[int]) -> List[int]:
        result = [0]*len(nums)
        i,j,k=0,len(nums)-1,len(nums)-1
        while i<=j:
            if nums[i] ** 2 < nums[j]**2:
                result[k] = nums[j]**2
                j-=1
            else:
                result[k] = nums[i]**2
                i+=1
            k-=1
        return result
            
class TestSolution(unittest.TestCase):
    def setUp(self):
        self.solution = Solution()

    def test_sortedSquares(self):
        # 测试用例1：包含负数、零和正数
        self.assertEqual(self.solution.sortedSquares([-4, -1, 0, 3, 10]), [0, 1, 9, 16, 100])
        # 测试用例2：全为负数
        self.assertEqual(self.solution.sortedSquares([-7, -3, -1]), [1, 9, 49])
        # 测试用例3：全为正数
        self.assertEqual(self.solution.sortedSquares([2, 3, 5]), [4, 9, 25])
        # 测试用例4：包含重复元素
        self.assertEqual(self.solution.sortedSquares([-2, -2, 2, 2]), [4, 4, 4, 4])
        # 测试用例5：空列表
        self.assertEqual(self.solution.sortedSquares([]), [])

if __name__ == '__main__':
    unittest.main()