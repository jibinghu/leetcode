#include "../definition.cpp"
#include <cmath>

class Solution{
    public:
        ListNode* getIntersectionNode(ListNode* headA, ListNode* headB){
            // 注意题目描述中的关键点：在intersection之后的链表都是一致的
            ListNode* currentA = headA;
            ListNode* currentB = headB;
            int numA=0,numB=0;
            while(currentA != nullptr){
                currentA = currentA->next;
                numA++;
            }
            while(currentB != nullptr){
                currentB = currentB->next;
                numB++;
            }
            currentA = headA;
            currentB = headB;
            int recursion = abs(numA - numB);
            while(recursion--){
                currentA = currentA->next;
            }
            while(currentA != nullptr && currentB != nullptr){
                if(currentA==currentB)
                    return currentA;
                currentA = currentA->next;
                currentB = currentB->next;
            }
            return nullptr;
        }
};