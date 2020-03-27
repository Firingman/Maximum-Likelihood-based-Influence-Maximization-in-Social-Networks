
import networkx as nx
import numpy as np
import random

datasets = []
f = open("1.txt")
data = f.read()
rows = data.split('\n')


for row in rows:
    split_row = row.split('\t')
    name = (int(split_row[0]), int(split_row[1]))
    datasets.append(name)
G = nx.DiGraph()
G.add_edges_from(datasets)


allnodes = G.nodes()
all_nodes = list(allnodes)
alledges = list(G.edges())
print(alledges)
print(all_nodes)

length = len(all_nodes)

T_values = [] # 每条边的权重的存储
# T = np.zeros(shape=(length,length))
#
# for i in range(length):
#     for j in range(length):
#         if (all_nodes[i],all_nodes[j]) in alledges:
#             # print(all_nodes[i], all_nodes[j])
#             T[i][j] = random.uniform(0, 0.09)
#         else:
#             T[i][j] = 0
#
# print('初始化的传播矩阵：') #这边不应该是对称的，有向网络
# print(T)
print(len(alledges))
for i in range(len(alledges)):
    a = random.uniform(0, 0.09)
    T_values.append(a)
print("每条边的权重值得列表：")
print(T_values)

# print("这些权重中的最小值：")
# print(min(T_values))
# print("这些权重中最小值的边：")
# print(T_values.index(min(T_values)))
# print("最小的值所在的边：")
# a = list(alledges[T_values.index(min(T_values))])[0]
# b = list(alledges[T_values.index(min(T_values))])[1]
#
# print(a)
# print(b)
#
# print(type(alledges[T_values.index(min(T_values))]))
# #初始化完成


'''
TCBS
'''
U = []#这个是存储0入度的节点

E_edge = [] #初始化E' = 空串

while len(all_nodes) > 0: #循环的条件
    node_zero_degree = []
    out_degrees = dict((u, 0) for u in G)
    in_degrees = dict((u, 0) for u in G)

    for u in G:
        out_degrees[u] = len(G[u])

    for u in G:
        for v in G[u]:
            in_degrees[v] += 1
    print("图中所有节点的入度：")
    print(in_degrees)
    print("图中所有节点的出度：")
    print(out_degrees)


    for i in in_degrees:
        # print("this_node:")
        # print(i)
        # print("the_in_degree_of_this_node:")
        # print(in_degrees[i])
        if in_degrees[i] == 0:
            node_zero_degree.append(i)
            print(i)
            all_nodes.remove(i)
        else:
            E_edge.append(alledges[T_values.index(min(T_values))]) #E加上需要去除掉的边
            T_values.remove(min(T_values)) #删除对应边上的权重
            alledges.remove(alledges[T_values.index(min(T_values))]) #删除原来所有上的边
            a = list(alledges[T_values.index(min(T_values))])[0]
            b = list(alledges[T_values.index(min(T_values))])[1]
            G.add_edges_from(alledges)


    print("节点入度为0：")
    print(node_zero_degree)
    U.append(node_zero_degree)
    print("这个是存储0入度的节点:")
    print(U)
    for i in node_zero_degree:
        l = list(G.neighbors(i))

        print(l)
        for k in l:
            print(k)
            E_edge.append((i, k))
            G.remove_edge(i, k) #删除跟这个入度为0的节点的相关的边 出度边

    print(E_edge)
    print(G.edges())
    print(list(G.edges()))
    for i in list(G.edges()):
        print(i)
        c = list(i)[0]
        d = list(i)[1]
        G.remove_edge(c, d)
    print(all_nodes)

    for i in node_zero_degree: #删除掉这个入度为0的节点
        G.remove_node(i)
    print(list(G.nodes))

print("节点入度为0的U：")
print(U)
print("所有去掉的边：")
print(E_edge)

L = []
for i in range(len(U)-1):
    for j in range(i+1, len(U)):
        sum1 = 0
        for k in U[i]:
            F = []
            for l in U[j]:
                if (k, l) in alledges:
                    sum1 = sum1 + (1 - T_values[alledges.index(k, l)])
        L.append(1 - sum1)

K_Location = L.index(min(L))
print("选择的节点：")
print(K_Location)








