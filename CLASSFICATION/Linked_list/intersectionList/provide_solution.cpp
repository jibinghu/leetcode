#include "../definition.cpp"
#include <cmath>

class Solution {
public:
    ListNode* getIntersectionNode(ListNode* headA, ListNode* headB) {
        // 检查空链表情况
        if (headA == nullptr || headB == nullptr) {
            return nullptr;
        }

        // 计算链表长度
        int lenA = getLength(headA);
        int lenB = getLength(headB);

        // 对齐链表
        ListNode* currentA = headA;
        ListNode* currentB = headB;
        int lengthDiff = abs(lenA - lenB);

        if (lenA > lenB) {
            currentA = moveForward(currentA, lengthDiff);
        } else {
            currentB = moveForward(currentB, lengthDiff);
        }

        // 查找交点
        while (currentA != nullptr && currentB != nullptr) {
            if (currentA == currentB) {
                return currentA;
            }
            currentA = currentA->next;
            currentB = currentB->next;
        }

        return nullptr; // 没有交点
    }

private:
    // 计算链表长度
    int getLength(ListNode* head) {
        int length = 0;
        while (head != nullptr) {
            length++;
            head = head->next;
        }
        return length;
    }

    // 将链表指针向前移动 steps 步
    ListNode* moveForward(ListNode* node, int steps) {
        while (steps-- > 0 && node != nullptr) {
            node = node->next;
        }
        return node;
    }
};