#include <iostream>
#include "definition.cpp"
using namespace std;

class Solution{
    public:
    ListNode* removeElements(ListNode* head, int val){
        // 添加虚拟头结点实现对头结点的删除
        ListNode* dummy = new ListNode(0);
        dummy->next = head;
        ListNode* cur = dummy;
        while (cur->next != nullptr){
            
            if(cur->next->val == val){
                cur->next = cur->next->next;
            }
            else{
                cur = cur->next;
            }
        }
        return dummy->next;
        delete dummy;
    }
};