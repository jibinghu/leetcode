#include <vector>
#include <iostream>
using namespace std;

class Solution
{
public:
    vector<vector<int>> generateMatrix(int n)
    {
        vector<vector<int>> result(n, vector<int>(n, 0));
        int top = 0, bottom = n - 1, left = 0, right = n - 1, count = 1;
        while (count <= n * n)
        {
            for (int i = left; i <= right; i++)
            {
                result[top][i] = count++;
            }
            top++;
            for (int j = top; j <= bottom; j++)
            {
                result[j][right] = count++;
            }
            right--;
            for (int k = right; k >= left; k--)
            {
                result[bottom][k] = count++;
            }
            bottom--;
            for (int m = bottom; m >= top; m--)
            {
                result[m][left] = count++;
            }
            left++;
        }
        return result;
    }
};

int main()
{
    int n;
    cin >> n;
    Solution solu;
    vector<vector<int>> vec = solu.generateMatrix(n);
    for (auto &i : vec)
    {
        for (auto &j : i)
        {
            cout << j << ' ';
        }
        cout << endl;
    }
}