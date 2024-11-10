#include <vector>
#include <iostream>
using namespace std;

class Solution {
public:
    int removeElement(vector<int>& nums, int val) {
        int left = 0, right = nums.size();
        while (left < right) {
            if (nums[left] == val) {
                for (int i = left + 1; i < right; i++) {
                    // 注意right=size，要避免数组的越界
                    nums[i - 1] = nums[i];
                }
                right--;
            } else {
                left++;
            }
        }
        return right;
    }
};


int main(int argc,char **argv){
    vector<int> vec={3,2,2,3,4,3,4,4,5};
    int val = 3;
    int con = removeElement(vec,val);
    for(int i=0;i<con;i++){
        printf("%d ",vec[i]);
    }
}