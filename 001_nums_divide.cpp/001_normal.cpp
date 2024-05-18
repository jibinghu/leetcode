#include <limits.h>
#include <stdexcept>

class Solution{
    public:
        int divide(int dividend, int divisor){

            // 异常处理
            if (divisor == 0)
                throw std::invalid_argument("Divisor can't be zero.\n");

            // 定义最值
            const int MAX_INT = INT_MAX;
            const int MIN_INT = INT_MIN;

            if (dividend == MIN_INT) {
                if (divisor == 1) {
                    return INT_MIN;
                }
                if (divisor == -1) {
                    return INT_MAX;
                }
        }

            // 溢出判断
            if (dividend == MAX_INT && divisor == -1)
                return MAX_INT;

            // 符号判断
            bool negative = (dividend < 0) != (divisor < 0);

            // 抽象为绝对值
            long long abs_dividend = std::abs(static_cast<long long>(dividend));
            long long abs_divisor = std::abs(static_cast<long long>(divisor));

            // 除法逻辑
            int result = 0;
            while(abs_dividend >= abs_divisor){
                long long tmp = abs_divisor, mutiple = 1;
                while(abs_dividend >= (tmp << 1)){
                    tmp <<= 1;
                    mutiple <<= 1;
                }
                abs_dividend -= tmp;
                result += mutiple; 
            }

            result = negative ? (-result) : result;
            return result;
        }
};

int main()
{
    Solution s1 = Solution();
    printf("%d",s1.divide(-2147483648,1));
}