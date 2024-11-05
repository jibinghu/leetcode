// O(n+nlogn)
#include <vector> 
#include <algorithm> 
#include <iostream> 
using namespace std; 

class Solution{
    public:
        vector<int> sortedSquares(vector<int>& A){
            for (int i=0;i<A.size();i++){
                A.at(i) *= A.at(i);
            }
            sort(A.begin(),A.end());
            return A;
        }
};

int main(int argc, char** argv){
    Solution solu;
    vector<int> nums = {-4,-1,0,3,10};
    vector<int> result = solu.sortedSquares(nums);
    for(int res:result)
        cout << res << ' ';
    cout << endl;
}