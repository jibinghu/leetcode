#include <iostream>
#include <stack>
#include <string>
#include <cctype>

// 定义二叉树节点结构
struct TreeNode {
    char value; // 操作符或操作数
    TreeNode* left;
    TreeNode* right;
    TreeNode(char val) : value(val), left(nullptr), right(nullptr) {}
};

// 判断字符是否为操作符
bool isOperator(char c) {
    return c == '+' || c == '-' || c == '*' || c == '/';
}

// 获取操作符的优先级
int precedence(char op) {
    if (op == '+' || op == '-') return 1;
    if (op == '*' || op == '/') return 2;
    return 0;
}

// 将中缀表达式转换为后缀表达式
std::string infixToPostfix(const std::string& infix) {
    /*
    中缀表达式转后缀表达式关键点：
    - 需要一个栈和一个string表达式 -
    1. 遇到数字直接输出；
    2. 遇到左括号直接进栈；
    3. 遇到右括号将栈内左括号之前的所有操作符弹出；
    4. 遇到操作符，将栈内所有优先级高的操作符弹出后再进栈。
    - 栈的方法：
    1. stack.top();
    2. stack.pop();
    3. stack.push();
    4. stakc.empty();
    5. stack.
    */ 
    std::stack<char> operators;
    std::string postfix;
    for (char ch : infix) {
        if (std::isdigit(ch)) {
            postfix += ch;
        } else if (ch == '(') {
            operators.push(ch);
        } else if (ch == ')') {
            while (!operators.empty() && operators.top() != '(') {
                postfix += operators.top();
                operators.pop();
            }
            if (!operators.empty()) operators.pop(); // 弹出 '('
        } else if (isOperator(ch)) {
            while (!operators.empty() && precedence(operators.top()) >= precedence(ch)) {
                postfix += operators.top();
                operators.pop();
            }
            operators.push(ch);
        }
    }
    while (!operators.empty()) {
        postfix += operators.top();
        operators.pop();
    }
    return postfix;
}

// 根据后缀表达式构建表达式二叉树
TreeNode* constructTree(const std::string& postfix, int& index) {
    /*
    构建表达式二叉树一般都是用递归的方法(从后往前递归构建二叉树)：
    如果当前字符是操作符，创建节点 node，然后递归构建 right 和 left 子节点，表示该操作符的操作数；
	如果当前字符是操作数，直接返回对应的节点
    */
    if (index < 0) return nullptr;
    char ch = postfix[index--];
    TreeNode* node = new TreeNode(ch);
    if (isOperator(ch)) {
        node->right = constructTree(postfix, index);
        node->left = constructTree(postfix, index);
    }
    return node;
}

// 构建表达式二叉树的接口函数
TreeNode* buildExpressionTree(const std::string& infix) {
    std::string postfix = infixToPostfix(infix);
    std::cout << "对应的后缀表达式：" << postfix << std::endl;
    int index = postfix.size() - 1;
    return constructTree(postfix, index);
}

// 中序遍历表达式二叉树
void inorderTraversal(TreeNode* root) {
    if (root) {
        if (root->left) {
            std::cout << "(";
            inorderTraversal(root->left);
        }
        std::cout << root->value;
        if (root->right) {
            inorderTraversal(root->right);
            std::cout << ")";
        }
    }
}

// 释放二叉树内存
void freeTree(TreeNode* root) {
    if (root) {
        freeTree(root->left);
        freeTree(root->right);
        delete root;
    }
}

int main() {
    std::string infixExpression;
    std::cout << "请输入中缀表达式（仅包含数字和操作符 +, -, *, /，不含空格）：";
    std::cin >> infixExpression;

    TreeNode* expressionTree = buildExpressionTree(infixExpression);

    std::cout << "中序遍历结果：";
    inorderTraversal(expressionTree);
    std::cout << std::endl;

    freeTree(expressionTree);
    return 0;
}
