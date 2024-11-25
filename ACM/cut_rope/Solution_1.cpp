#include <iostream>
#include <iomanip> // 用于设置输出精度
#include <vector>
#include <algorithm> // 用于 std::max 和 std::min

using namespace std;

// 函数用于判断是否能通过当前长度maxLength切割出K条绳子
bool canCut(const vector<double>& ropes, int K, double maxLength) {
    int count = 0; // 记录切割出的绳子条数
    for (double rope : ropes) {
        count += static_cast<int>(rope / maxLength); // 每根绳子可以切出的条数
    }
    return count >= K; // 判断是否至少切出K条
}

double findMaxLength(const vector<double>& ropes, int K) {
    double low = 0.0; // 最小可能长度
    double high = *max_element(ropes.begin(), ropes.end()); // 所有绳子中最长的那一条
    double result = 0.0; // 保存结果

    // 二分查找精确到0.01
    while (high - low > 1e-3) { // 精度控制
        double mid = (low + high) / 2; // 当前中间长度
        if (canCut(ropes, K, mid)) {
            result = mid; // 更新结果
            low = mid;    // 尝试更长的绳子
        } else {
            high = mid;   // 尝试更短的绳子
        }
    }
    return result;
}

int main() {
    int N, K;
    // 输入绳子数量和需要切割的绳子数量
    cin >> N >> K;
    vector<double> ropes(N);

    // 输入每条绳子的长度
    for (int i = 0; i < N; ++i) {
        cin >> ropes[i];
    }

    // 调用函数求解最大长度
    double maxLength = findMaxLength(ropes, K);

    // 设置输出精度为两位小数
    cout << fixed << setprecision(2) << maxLength << endl;

    return 0;
}