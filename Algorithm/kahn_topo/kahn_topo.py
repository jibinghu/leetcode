from collections import deque

def kahn_topological_sort(graph):
    # 计算每个节点的入度
    in_degree = {node: 0 for node in graph}

    for node in graph:
        for neighbor in graph[node]:
            in_degree[neighbor] += 1

    # 初始化队列，将入度为 0 的节点加入队列
    queue = deque([node for node in graph if in_degree[node] == 0])
    topological_order = []

    # 进行拓扑排序
    while queue:
        node = queue.popleft()
        topological_order.append(node)
        
        # 遍历该节点的邻接节点，并更新入度
        for neighbor in graph[node]:
            in_degree[neighbor] -= 1
            # 如果邻接节点的入度变为 0，则加入队列
            if in_degree[neighbor] == 0:
                queue.append(neighbor)
    
    # 检查是否所有节点都在拓扑排序中（检测环）
    if len(topological_order) == len(graph):
        return topological_order
    else:
        return "Graph has a cycle, topological sorting is not possible."

# 示例图（邻接表表示）
graph = {
    "起床": ["洗脸", "刷牙"],
    "洗脸": ["早餐", "刷牙"],
    "刷牙": ["早餐"],
    "早餐": []
}

print(kahn_topological_sort(graph))