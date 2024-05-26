
#include <vetor>
class Solution {
public:
    int trap(vector<int>& height) {
        // 按行进行计算
        int n = height.size();
        int capacity = 0;
        stack<int> stk;

        for (int i=0; i<n; i++){
            while(stk.size() && height[i] >= height[stk.top()]){
                int l = height[stk.top()];
                stk.pop();
                if (stk.size()){
                    int h = min(height[i], height[stk.top()]);
                    capacity += (h-l) * (i - stk.top() - 1);
                }
            }
            stk.push(i);
        }
        return capacity;
    }
};