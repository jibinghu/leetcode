#include <vector>

class Solution {
public:
    std::vector<int> countBits(int n) {
        std::vector<int> ans;
        for(int i=0;i<n+1;i++){
            int count = __builtin_popcount(i);
            ans.push_back(count);
        }
        return ans;
    }
};