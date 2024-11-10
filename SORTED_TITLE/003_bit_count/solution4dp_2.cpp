#include <vector>

class Solution {
public:
    std::vector<int> countBits(int n) {
        // 让我们从更直观的最低有效位来看
        std::vector<int> bits(n + 1);
        for (int i = 1; i <= n; i++) {
            // 每个数是与>>1的数进行比较
            bits[i] = bits[i >> 1] + (i & 1);
        }
        return bits;
    }
};

// class Solution {
// public:
//     std::vector<int> countBits(int n) {
//         std::vector<int> ans;
//         if (n == 0) {
//             return {0};
//         }
//         ans = {0, 1};
//         for (int i = 2; i <= n; ++i) {
//             int t1 = i / 2;
//             int t2 = i % 2;
//             ans.push_back(ans[t1] + t2);
//         }
//         return ans;
//     }
// };