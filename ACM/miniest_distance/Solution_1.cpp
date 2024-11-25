#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

int minDistance(string word1, string word2) {
    int m = word1.size();
    int n = word2.size();

    // 创建 DP 表格，大小为 (m+1) x (n+1)
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

    // 初始化第一行和第一列
    for (int i = 0; i <= m; ++i) {
        dp[i][0] = i; // 删除操作
    }
    for (int j = 0; j <= n; ++j) {
        dp[0][j] = j; // 插入操作
    }

    // 填充 DP 表格
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (word1[i - 1] == word2[j - 1]) {
                // 当前字符相等，不需要额外操作
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                // 三种操作的最小值 + 1
                dp[i][j] = min({dp[i - 1][j],    // 删除操作
                                dp[i][j - 1],    // 插入操作
                                dp[i - 1][j - 1] // 替换操作
                               }) + 1;
            }
        }
    }

    // 返回最终结果
    return dp[m][n];
}

int main() {
    string word1, word2;
    cout << "Enter the first word: ";
    cin >> word1;
    cout << "Enter the second word: ";
    cin >> word2;

    int result = minDistance(word1, word2);
    cout << "The minimum edit distance is: " << result << endl;

    return 0;
}