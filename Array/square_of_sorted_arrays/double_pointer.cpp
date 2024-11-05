#include <vector>
#include <iostream>
using namespace std;
// 双指针遍历 O(n)
class Solution{
    public:
        vector<int> sortedSquares(vector<int>& A){
            vector<int> result(A.size());
            int i=0,j=A.size()-1,k=A.size()-1;
            while(i<=j && k>=0){ // 注意 i<=j，对于最后一个元素的处理
                if(A.at(i)*A.at(i)<A.at(j)*A.at(j)){
                    result[k--] = A.at(j)*A.at(j);
                    j--;
                }
                else{
                    result[k--] = A.at(i)*A.at(i);
                    i++;
                }
                // k--;
            }
            return result;
        }
};

int main(int argc, char** argv){
    Solution solu;
    vector<int> nums = {-4,-1,0,3,10};
    vector<int> result = solu.sortedSquares(nums);
    for(int res:result)
        cout << res << ' ';
}