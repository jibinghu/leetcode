#include <vector>
#include <climits>
#include <iostream>
using namespace std;
// O(n)，两次循环，只不过第二次循环是在滑动窗口里，记为O(2*n);
class Solution{
    public:
        int minSubArrayLen(int s, vector<int>& nums){
            int i=0,sum=0,current=0,result=INT_MAX;
            for (int j=0; j<nums.size();j++){
                sum += nums.at(j);
                while(sum >= s){
                    // 注意这里的i并没有每次回溯到0位置
                    current = j - i + 1;
                    result = result > current ? current : result;
                    sum -= nums.at(i++);
                }
            }
            return result == INT_MAX ? 0 : result;
        }
};

int main(){
    int s = 7;
    vector<int> nums={2,3,1,2,4,3};
    Solution solu;
    cout << solu.minSubArrayLen(s,nums);
    cout << endl;
}