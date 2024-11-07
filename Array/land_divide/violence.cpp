#include <vector>
#include <iostream>
#include <climits>
using namespace std;
// O(n^3)
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
    // 逻辑
    // row-split
    int sign=INT_MAX;
    for(int s=0;s<n-1;s++){
        int asum=0,bsum=0;
        for(int i=0;i<=s;i++){
            for(int j=0;j<m;j++){
                asum+=land[i][j];
            }
        }
        for(int i=s+1;i<=n-1;i++){
            for(int j=0;j<m;j++){
                bsum+=land[i][j];
            }
        }
        sign = abs(asum - bsum) > sign ? sign : abs(asum - bsum);
    }
    // column-split
    for(int s=0;s<m-1;s++){
        int asum=0,bsum=0;
        for(int i=0;i<n;i++){
            for(int j=0;j<=s;j++){
                asum+=land[i][j];
            }
        }
        for(int i=0;i<n;i++){
            for(int j=s+1;j<m;j++){
                bsum+=land[i][j];
            }
        }
        sign = abs(asum - bsum) > sign ? sign : abs(asum - bsum);
    }
    cout << sign << endl;
    return 0;
}