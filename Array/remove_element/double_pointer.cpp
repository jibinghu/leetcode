#include <vector>
using namespace std;

class Solution {
public:
    int removeElement(vector<int>& nums, int val) {
        int left = 0; // 用于记录新数组的位置

        // 遍历数组
        for (int right = 0; right < nums.size(); right++) {
            // 如果当前元素不等于val，将其移到left指针位置
            if (nums[right] != val) {
                nums[left++] = nums[right];
            }
        }

        // 返回不包含指定值的新数组的长度
        return left;
    }
};