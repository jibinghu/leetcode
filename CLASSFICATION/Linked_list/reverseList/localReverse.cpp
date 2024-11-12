#include "../definition.cpp"

class Solution{
    public:
        // 原地翻转链表，避免内存浪费
        ListNode* reverseList(ListNode* head){
            if(head == nullptr || head->next == nullptr){
                return head;
            }
            else if(head->next->next == nullptr){
                ListNode* newNode = head->next;
                newNode->next = head;
                head->next = nullptr;
                return newNode;
            }
            else {
                ListNode* current = head->next;
                ListNode* previous = head;
                ListNode* next = head->next->next;

                while(next != nullptr){
                    current->next = previous;
                    previous = current;
                    current = next;
                    next = next->next;
                }
                current->next = previous;
                head->next = nullptr;
                return current;
            }

        }
};