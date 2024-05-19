#include <stdio.h>
#include <cmath>


void long_division(int dividend, int divisor, int *quotient, int *remainder) {
    *quotient = 0;
    *remainder = 0;

    int sign = 1;
    if ((dividend < 0 && divisor > 0) || (dividend > 0 && divisor < 0)) {
        sign = -1;
    }

    dividend = std::abs(dividend);
    divisor = abs(divisor);

    for (int i = 31; i >= 0; i--) {
        if ((*remainder << 1) + ((dividend >> i) & 1) >= divisor) {
            *remainder = (*remainder << 1) + ((dividend >> i) & 1) - divisor;
            *quotient = (*quotient << 1) | 1;
        } else {
            *remainder = (*remainder << 1) + ((dividend >> i) & 1);
            *quotient = *quotient << 1;
        }
    }

    *quotient *= sign;
    if (sign < 0) {
        *remainder = -(*remainder);
    }
}

int main() {
    int dividend = -1010369383;
    int divisor = -2147483648;
    int quotient, remainder;

    long_division(dividend, divisor, &quotient, &remainder);

    printf("Quotient: %d\n", quotient);
    printf("Remainder: %d\n", remainder);

    return 0;
}
