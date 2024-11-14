#include "../definition.cpp"

class Solution{
    public:
        ListNode* reverseList(ListNode* head){
            ListNode* previous = nullptr;
            ListNode* current = head;
            ListNode* tmp;

            while(current != nullptr){
                tmp = current->next;
                current->next = previous;
                previous = current;
                current = tmp;
            }
            return previous;
        }
};