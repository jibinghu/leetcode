#include <vector>

class Solution {
public:
    int counts(int x){
        int nums = 0;
        while (x > 0){
            x = x & (x - 1);
            nums++;
        }
        return nums;
    }

    std::vector<int> countBits(int n) {
        std::vector<int> vec;
        for (int i=0; i<n+1; i++){
            int k = counts(i);
            vec.push_back(k);
        }
        return vec;
    }
};