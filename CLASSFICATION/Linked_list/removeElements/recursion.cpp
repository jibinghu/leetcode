#include "../definition.cpp"
using namespace std;

class Solution {
public:
    ListNode* removeElements(ListNode* head, int val) {
        if (head == nullptr) {
            return nullptr; // 基础情况：链表为空，返回 nullptr
        }

        // 递归调用处理链表的下一个节点
        if (head->val == val) {
            ListNode* newNode = removeElements(head->next, val);
            delete head; // 删除当前节点以释放内存
            return newNode; // 返回新链表的头节点
        } else {
            head->next = removeElements(head->next, val);
            return head; // 返回当前节点作为新链表的头节点
        }
    }
};