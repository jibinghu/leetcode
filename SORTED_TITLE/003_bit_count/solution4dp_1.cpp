#include <vector>

class Solution
{
public:
    std::vector<int> countBits(int n)
    {
        // dp思想
        std::vector<int> bits(n + 1);
        bits[0] = 0;
        int highbit = 0;
        for (int i = 1; i <= n; i++)
        {
            // 只有2的m次幂的KB式子才会等于零，也就是二进制1的个数为一
            if (!(i & (i - 1)))
                highbit = i;
            bits[i] = bits[i - highbit] + 1;
        }
        return bits;
    }
};