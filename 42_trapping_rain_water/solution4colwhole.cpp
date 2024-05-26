#include <vector>

class Solution {
public:
    int trap(vector<int>& height) {
        // 双指针按列计算
        int n = height.size();
        int left = 0, right = n-1;
        int leftMax = 0, rightMax = 0;
        int capacity=0;
        while(left < right){
            leftMax = max(height[left], leftMax);
            rightMax = max(height[right], rightMax);
            if (leftMax < rightMax){
               capacity += leftMax - height[left]; 
               left++;
            }
            else{
                capacity += rightMax - height[right];
                right--;
            }
        }
        return capacity;
    }
};