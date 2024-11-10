#include <vector>
#include <iostream>
#include <climits>
using namespace std;
// O(n^2)
int main(){
    // 输入
    int n,m;
    cin >> n >> m;
    vector<vector<int>> land(n, vector<int>(m,0));
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            cin >> land[i][j];
        }
    } 
    // 行/列前缀和
    vector<int> row_pre(n);
    vector<int> col_pre(m);
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            row_pre[i] += land[i][j];
        }
    }
    for(int j=0;j<m;j++){
        for(int i=0;i<n;i++){
            col_pre[j] += land[i][j];
        }
    }
    // 逻辑
    // row-split
    int sign=INT_MAX;
    for(int s=0;s<n-1;s++){
        int asum=0,bsum=0;
        for(int i=0;i<=s;i++){
            asum+=row_pre[i];
        }
        for(int i=s+1;i<=n-1;i++){
            bsum+=row_pre[i];
        }
        sign = abs(asum - bsum) > sign ? sign : abs(asum - bsum);
    }
    // column-split
    for(int s=0;s<m-1;s++){
        int asum=0,bsum=0;
        for(int j=0;j<=s;j++){
            asum+=col_pre[j];
        }
        for(int j=s+1;j<m;j++){
            bsum+=col_pre[j];
        }
        sign = abs(asum - bsum) > sign ? sign : abs(asum - bsum);
    }
    cout << sign << endl;
    return 0;
}