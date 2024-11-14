#include "../definition.cpp"

/*
 * 对链表操作时，尽量直接设置一个虚拟头结点
 * 注意判断和合并循环边界条件
 * next 可以指向空但不能为空
 */

class Solution{
    public:
        ListNode* swapPairs(ListNode* head){
            // 设置一个虚拟头结点
            ListNode* dummyNode = new ListNode(0);
            dummyNode->next = head;
            // 分三步解决问题
            /* 1，2，3，4
             * 第一：1 -> 3
             * 第二：3 -> 2
             * 第三：2 -> 4
            */
           ListNode* current = dummyNode;
           while(current->next != nullptr && current->next->next != nullptr){
                ListNode* ltmp = current->next;
                ListNode* rtmp = ltmp->next->next;

                current->next = ltmp->next; // step 1
                current->next->next = ltmp; // step 2
                ltmp->next = rtmp; // step 3

                current = current->next->next; // 后边已经排好了，直接向后移动两位
           }
           return dummyNode->next;
        }
};