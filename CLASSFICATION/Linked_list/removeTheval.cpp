#include "definition.cpp"
using namespace std;

class Solution{
    public:
    ListNode* removeElements(ListNode* head, int val){
        // 对于头结点的删除单独设置逻辑
        // 这里对于每次删除的结点，需要设置tmp临时结点清理内存
        while (head != nullptr && head->val == val){
            head = head->next;
        }
        ListNode* cur = head;
        while (cur != nullptr && cur->next != nullptr){
            if(cur->next->val == val){
                cur->next = cur->next->next;
            }
            else{
                cur = cur->next;
            }
        }
        return head;
    }
};