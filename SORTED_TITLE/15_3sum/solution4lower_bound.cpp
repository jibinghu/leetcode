#include <vector>
using namespace std;

class Solution {
public:
    vector<vector<int>> threeSum(vector<int>& nums) {
        vector<vector<int>> res;
        int n = nums.size();
        sort(nums.begin(), nums.end());
        for (int i = 0; i < n - 2; i++) {
            if (i > 0 && nums[i] == nums[i - 1])
                continue;
            for (int j = i + 1; j < n - 1; j++) {
                if (j > i + 1 && nums[j] == nums[j - 1])
                    continue;
                int endnum = 0 - nums[i] - nums[j];
                auto it = lower_bound(nums.begin(), nums.end(), endnum);
                if (it != nums.end() && *it == endnum){
                    res.push_back({nums[i],nums[j],endnum});
                }
            }
        }
        return res;
    }
};