class Solution {
public:
    vector<vector<int>> threeSum(vector<int>& nums) {
        vector<vector<int>> res;
        int n = nums.size();
        sort(nums.begin(), nums.end());
        for (int i = 0; i < n - 2; i++) {
            if (i > 0 && nums[i] == nums[i - 1])
                continue;
            if (nums[i+1]+nums[i+2]>-nums[i]){
                return res;
            if (nums[n-1]+nums[n]<-nums[i]){
                continue;
            }
            }
            int left=i+1,right=n-1;
            while(left < right){
                if (nums[left]+nums[right]==-nums[i]){
                    res.push_back({nums[i],nums[left],nums[right]});
                    left++;
                    while(left<right && nums[left]==nums[left-1]) left++;
                    right--;
                    while(left<right && nums[right]==nums[right+1]) right--;
                }
                else if (nums[left]+nums[right]<-nums[i])
                    left++;
                else
                    right--;
            }
        }
        return res;
    }
};