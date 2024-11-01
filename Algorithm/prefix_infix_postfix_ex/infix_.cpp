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
