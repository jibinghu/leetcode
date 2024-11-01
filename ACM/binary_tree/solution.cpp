#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// 构建二叉树节点
struct TreeNode {
    int val;
    TreeNode* left;
    TreeNode* right;
    TreeNode(int x) : val(x), left(NULL), right(NULL) {}
};

// 建树函数，根据数组建立二叉树
TreeNode* buildTree(const vector<int>& nodes, int index) {
    if (index >= nodes.size() || nodes[index] == 0) {
        return nullptr;
    }
    TreeNode* root = new TreeNode(nodes[index]);
    root->left = buildTree(nodes, 2 * index + 1);
    root->right = buildTree(nodes, 2 * index + 2);
    return root;
}

// 递归函数计算最少基站数量
int minBaseStations(TreeNode* root, int& count) {
    if (!root) return 2; // 空节点返回2，表示不需要基站
    int left = minBaseStations(root->left, count);
    int right = minBaseStations(root->right, count);
    if (left == 0 || right == 0) {
        count++;
        return 1; // 当前节点需要建设基站
    }
    return (left == 1 || right == 1) ? 2 : 0; // 如果子节点有基站，当前节点不需要基站
}

int main() {
    string input;
    cin >> input;
    vector<int> nodes(input.size());
    for (int i = 0; i < input.size(); ++i) {
        nodes[i] = input[i] - '0';
    }

    TreeNode* root = buildTree(nodes, 0);
    int count = 0;
    if (minBaseStations(root, count) == 0) {
        count++;
    }

    cout << count << endl;
    return 0;
}