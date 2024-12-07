#include <string>
using namespace std;

class Solution {
public:
    bool isAnagram(string s, string t) {
        // 解法确实牛，这个字母映射到数字的技巧可以经常用到
        // 另外数组在C++中也是unordered_map的特例，搜索和增删改效率是O(1)
        int record[26] = {0};
        for(int i=0;i<s.size();i++){
            record[s[i]-'a']++;
        }
        for(int j=0;j<t.size();j++){
            record[t[j]-'a']--;
        }
        for(int k=0;k<26;k++){
            if(record[k]!=0)
                return false;
        }
        return true;
    }
};