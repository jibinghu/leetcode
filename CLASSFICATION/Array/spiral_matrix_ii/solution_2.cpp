#include <vector>
#include <iostream>
using namespace std;

class Solution {
public:
    vector<vector<int>> generateMatrix(int n) {
        vector<vector<int>> result(n, vector<int>(n, 0));
        int count = 1;  // 要填入的数字
        int start = 0;  // 起始点
        int size = n;   // 当前矩阵的大小

        while (size > 0) {
            // 上边从左到右填充
            for (int i = start; i < start + size; i++) {
                result[start][i] = count++;
            }
            // 右边从上到下填充
            for (int i = start + 1; i < start + size; i++) {
                result[i][start + size - 1] = count++;
            }
            // 下边从右到左填充（仅当 size > 1 时）
            if (size > 1) {
                for (int i = start + size - 2; i >= start; i--) {
                    result[start + size - 1][i] = count++;
                }
            }
            // 左边从下到上填充（仅当 size > 1 时）
            if (size > 1) {
                for (int i = start + size - 2; i > start; i--) {
                    result[i][start] = count++;
                }
            }

            // 缩小边界
            start++;  // 向内收缩
            size -= 2;  // 减少矩阵大小
        }

        return result;
    }
};

int main(){
    int n;
    cin >> n;
    Solution solu;
    vector<vector<int>> vec = solu.generateMatrix(n);
    for(auto& i:vec){
        for(auto& j:i){
            cout << j << ' ';
        }
        cout << endl;
    }
}