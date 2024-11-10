#include <string>
#include <iostream>

using namespace std;

class Solution
{
public:
    int charToInt(char cc)
    {
        return cc - '0';
    }
    char intToChar(int ii)
    {
        return ii + 48;
    }
    string addBinary(string a, string b)
    {
        int a_length = a.length() - 1;
        int b_length = b.length() - 1;
        int carry = 0;
        string result = "";

        while (a_length >= 0 || b_length >= 0 || carry)
        {
            carry += (a_length >= 0 ? charToInt(a[a_length]) : 0) +
                     (b_length >= 0 ? charToInt(b[b_length]) : 0);
            int bitty = carry % 2;
            carry = carry / 2;

            string::iterator it = result.begin();
            result.insert(it, intToChar(bitty));
            a_length--;
            b_length--;
        }
        return result;
    }
};

int main()
{
    Solution s1;
    string a = "101";
    string b = "10";
    printf("%s\n",s1.addBinary(a, b).c_str());
    puts(s1.addBinary(a,b).c_str());
    std::cout << s1.addBinary(a,b);
}