#include <vector>

class Solution {
public:
    std::vector<int> countBits(int n) {
        std::vector<int> ans;
        for(int i=0;i<n+1;i++){
            // CPP的内置函数__builtin_popcount()用于计算给定的整数的二进制表示中的 111 的数目
            int count = __builtin_popcount(i);
            ans.push_back(count);
        }
        return ans;
    }
};