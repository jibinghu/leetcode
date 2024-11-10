class Solution {
public:
    int removeElement(vector<int>& nums, int val) {
        int left = 0, right = nums.size();
        for (; left < right; left++) {
            if (nums[left] == val) {
                for (int i = left + 1; i < right; i++) {
                    nums[i - 1] = nums[i];
                }
                right--;
                left--;
            }
        }
        return right;
    }
};
