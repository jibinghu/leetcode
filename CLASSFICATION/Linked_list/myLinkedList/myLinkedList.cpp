/*
 * get(index)：获取链表中第 index 个节点的值。如果索引无效，则返回-1。
 * addAtHead(val)：在链表的第一个元素之前添加一个值为 val 的节点。插入后，新节点将成为链表的第一个节点。
 * addAtTail(val)：将值为 val 的节点追加到链表的最后一个元素。
 * addAtIndex(index,val)：在链表中的第 index 个节点之前添加值为 val  的节点。如果 index 等于链表的长度，则该节点将附加到链表的末尾。如果 index 大于链表长度，则不会插入节点。如果index小于0，则在头部插入节点。
 * deleteAtIndex(index)：如果索引 index 有效，则删除链表中的第 index 个节点。
*/
class myLinkedList{
    public:
        // 链表类内嵌套一个结点 struct
        struct ListNode
        {
            int val;
            ListNode* next;
            ListNode(int x): val(x), next(nullptr){};
        };
        // 类初始化
        myLinkedList(){
            _dummyNode = new ListNode(0);
            _size = 0;
        }
        // 获取链表中第 index 个节点的值。如果索引无效，则返回-1
        int get(int index){
            if(index > _size-1 || index < 0){
                return -1;
            }
            ListNode* cur = _dummyNode->next;
            while(index>0){
                cur = cur->next;
                index--;
            }
            return cur->val;
        }
        // 在链表的第一个元素之前添加一个值为 val 的节点。插入后，新节点将成为链表的第一个节点
        void addAtHead(int val){
            ListNode* newNode = new ListNode(val);
            newNode->next = _dummyNode->next;
            _dummyNode->next = newNode;
            _size++;
        }
        // 将值为 val 的节点追加到链表的最后一个元素
        void addAtTail(int val){
            ListNode* newNode = new ListNode(val);
            ListNode* cur = _dummyNode;
            while(cur->next != nullptr){
                cur = cur->next;
            }
            cur->next = newNode;
            _size++;
        }
        // 在链表中的第 index 个节点之前添加值为 val  的节点。如果 index 等于链表的长度，
        // 则该节点将附加到链表的末尾。如果 index 大于链表长度，则不会插入节点。如果index小于0，则在头部插入节点
        void addAtIndex(int index,int val){
            if(index<=0){
                addAtHead(val);
            }
            else if(index == _size){
                addAtTail(val);
            }
            else if(index>0 && index<_size){
                ListNode* newNode = new ListNode(val);
                ListNode* cur = _dummyNode;
                // 移动到插入点的前一个节点
                // 这里index从0开始
                while(index--){
                    cur = cur->next;
                }
                newNode->next = cur->next;
                cur->next = newNode;
                _size++;
            }
        }
        // 如果索引 index 有效，则删除链表中的第 index 个节点
        void deleteAtIndex(int index){
            if(index<0 || index>=_size)
                return;
            ListNode* cur = _dummyNode;
            while(index--){
                cur = cur->next;
            }
            cur->next = cur->next->next;
            _size--;
        }
    private:
        // 定义虚拟头结点和链表长度
        ListNode* _dummyNode;
        int _size;
};