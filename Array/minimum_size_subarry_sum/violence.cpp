#include <vector>
#include <iostream>
#include <climits>
using namespace std;
// O(n^2)
class Solution{
    public:
        int minSubArryLen(int s, vector<int>& nums){
            int sign=INT_MAX,current=0,sum=0;
            for (int i=0;i<nums.size();i++){
                sum = nums.at(i);
                if(sum >= s)
                    return current;
                for (int j=i+1;j<nums.size();j++){
                    sum += nums.at(j);
                    if (sum>=s){
                        current = j - i + 1;
                        sign = sign > current ? current : sign;
                        break;
                    }
                }
            }
            return sign == INT_MAX ? 0 : sign;
        }   
};

int main(){
    int s = 7;
    vector<int> nums={2,3,1,2,4,3};
    Solution solu;
    cout << solu.minSubArryLen(s,nums);
    cout << endl;
}