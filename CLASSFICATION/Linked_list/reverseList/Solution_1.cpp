#include "../definition.cpp"


// 双指针的经典应用，如果要删除倒数第n个节点，让fast移动n步，然后让fast和slow同时移动，直到fast指向链表末尾。删掉slow所指向的节点就可以了
class Solution{
    public:
        ListNode* removeNthFromEnd(ListNode* head, int n){
            ListNode* dummyNode = new ListNode(0);
            dummyNode->next = head;
            ListNode* fast = dummyNode;
            ListNode* slow = dummyNode;

            // 很简单的道理，当时没有想出来😭
            while(n--){
                fast = fast->next;
            }
            // 当fast到最后一个结点是，slow也到了倒数第n个的前一个结点
            while(fast->next != nullptr){
                fast = fast->next;
                slow = slow->next;
            }
            slow->next = slow->next->next;
            return dummyNode->next;
        }
};