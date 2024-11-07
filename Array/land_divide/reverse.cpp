#include <vector>
#include <iostream>
#include <climits>
#include <algorithm>
using namespace std;
// O(n^2)
// 逆向思维
int main(){
    // 输入
    int n,m,sum=0;
    cin >> n >> m;
    vector<vector<int>> land(n, vector<int>(m,0));
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            cin >> land[i][j];
            // 在这里直接计算出累加和
            sum += land[i][j];
        }
    } 
    // 逻辑
    // row-split
    int sign=INT_MAX,count=0;
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            count += land[i][j];
            if(j == m-1)
                sign = min(sign, abs(sum - count - count));
        }
    }
    // column-split
    count = 0;
    for(int j=0;j<m;j++){
        for(int i=0;i<n;i++){
            count += land[i][j];
            if(i == n-1)
                sign = min(sign, abs(sum - count - count));
        }
    }
    cout << sign << endl;
    return 0;
}