#include <vector>
using namespace std;

class Solution {
public:
    vector<int> countBits(int n) {
        // 让我们从更直观的最低有效位来看
        vector<int> bits(n + 1);
        for (int i = 1; i <= n; i++) {
            // 将最高有效位清零再加一
            bits[i] = bits[i & (i - 1)] + 1;
        }
        return bits;
    }
};