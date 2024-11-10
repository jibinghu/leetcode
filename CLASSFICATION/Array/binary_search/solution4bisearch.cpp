#include <vector>
using namespace std;
// O(logN) O(1)
/*二分法的前提是数组为有序数组，同时题目还强调数组中无重复元素*/

class Solution {
public:
    int search(vector<int>& nums, int target) {
        // 二分查找
        int n =nums.size();
        int left=0,right=n-1;
        while (left <= right){
            int mid = left + (right - left) >> 1; // 避免溢出，等同于(left+right)/2
            if (nums[mid] == target){
                return mid;
            }
            else if (nums[mid] < target){
                left = mid+1;
            }
            else {
                right = mid-1;
            }
        }
        return -1;
    }
};