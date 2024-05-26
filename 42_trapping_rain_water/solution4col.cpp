#include <vector>

class Solution {
public:
    int trap(vector<int>& height) {
        // 每列存储的雨水是左右最高列中最小值与当前列值差
        int n = height.size();

        vector<int> leftMax(n);
        vector<int> rightMax(n);

        leftMax[0] = height[0];
        for (int i=1; i<n; i++){
            leftMax[i] = max(leftMax[i-1], height[i]);
        }
        rightMax[n-1] = height[n-1];
        for (int j=n-2; j>0; j--){
            rightMax[j] = max(rightMax[j+1], height[j]);
        }

        int capacity = 0;
        for (int k=1; k<n-1; k++){
            capacity += min(leftMax[k], rightMax[k]) - height[k];
        }
        return capacity;
    }
};