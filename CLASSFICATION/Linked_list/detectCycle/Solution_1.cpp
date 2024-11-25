#include "../definition.cpp"

class Solution{
    public:
    ListNode *detectCycle(ListNode* head){
        if (head == nullptr) return nullptr;
        // 判断有无环：定义快慢指针fast和slow，他们两个的相遇位置一定是在环内，
        // 且slow指针一定还没有走完一圈，fast指针走了n圈（fast相对slow每次走一格）
        ListNode* fast = head;
        ListNode* slow = head;
        while(fast->next!=nullptr && fast->next->next!=nullptr){
            fast = fast->next->next;
            slow = slow->next;
            if(fast == slow) break;
        }
        if(fast->next == nullptr || fast->next->next == nullptr)
            return nullptr;
        // inedex1从头结点开始走，index2从相遇结点开始走，他们俩的相遇点就是入口结点
        ListNode* index1 = head;
        ListNode* index2 = fast;
        while(index1!=index2){
            index1 = index1->next;
            index2 = index2->next;
        }
        return index1;
    }
};