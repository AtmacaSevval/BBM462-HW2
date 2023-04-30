import networkx as nx
from itertools import combinations
import numpy as np
import operator
import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import table


G1 = nx.read_weighted_edgelist('ENZYMES_g106.edges')
G2 = nx.read_weighted_edgelist('aves-barn-swallow-contact-network.edges')


def create_random_graphs(G):
    degrees = [d for n, d in G.degree()]

    random_graphs = []
    for i in range(100):
        random_graph = nx.configuration_model(degrees, seed=i)
        random_graph = nx.Graph(random_graph)
        random_graphs.append(random_graph)

    return random_graphs


randomG1 = create_random_graphs(G1)
randomG2 = create_random_graphs(G2)

M1 = nx.Graph([(0, 1), (0, 3), (0, 2)])
M2 = nx.Graph([(0, 1), (0, 2), (1, 3)])
M3 = nx.Graph([(0, 1), (0, 3), (0, 2), (2, 3)])
M4 = nx.Graph([(0, 1), (0, 2), (2, 3), (1, 3)])
M5 = nx.Graph([(0, 1), (0, 2), (2, 3), (1, 3), (0, 3)])
M6 = nx.Graph([(0, 1), (0, 2), (2, 3), (1, 3), (0, 3), (1, 2)])

M3_1 = nx.Graph([(0, 1), (0, 2)])
M3_2 = nx.Graph([(0, 1), (1, 2), (0, 2)])


def find_motifs_4(graph):
    counts = [0, 0, 0, 0, 0, 0]

    for nodes in combinations(graph.nodes(), 4):
        subgraph = graph.subgraph(nodes)
        if len(subgraph) > 0 and nx.is_connected(subgraph):
            if nx.is_isomorphic(subgraph, M1):
                counts[0] += 1
            elif nx.is_isomorphic(subgraph, M2):
                counts[1] += 1
            elif nx.is_isomorphic(subgraph, M3):
                counts[2] += 1
            elif nx.is_isomorphic(subgraph, M4):
                counts[3] += 1
            elif nx.is_isomorphic(subgraph, M5):
                counts[4] += 1
            elif nx.is_isomorphic(subgraph, M6):
                counts[5] += 1

    return counts


def find_motifs_3(graph):
    counts = [0, 0]

    for nodes in combinations(graph.nodes(), 3):
        subgraph = graph.subgraph(nodes)
        if len(subgraph) > 0 and nx.is_connected(subgraph):
            if nx.is_isomorphic(subgraph, M3_1):
                counts[0] += 1
            elif nx.is_isomorphic(subgraph, M3_2):
                counts[1] += 1
    return counts


def find_average(randomGraph, divide):
    average_num_random_motif3 = [0, 0]
    average_num_random_motif4 = [0, 0, 0, 0, 0, 0]

    if divide == 1:
        average_num_random_motif3 = find_motifs_3(randomGraph)
        average_num_random_motif4 = find_motifs_4(randomGraph)

    else:
        for g in randomGraph:
            average_num_random_motif3 = list(map(operator.add, average_num_random_motif3, find_motifs_3(g)))
            average_num_random_motif4 = list(map(operator.add, average_num_random_motif4, find_motifs_4(g)))

    motif3 = np.divide(average_num_random_motif3, divide)
    motif4 = np.divide(average_num_random_motif4, divide)

    return motif3, motif4


average_random1_motif3, average_random1_motif4 = find_average(randomG1, 100)
average_real1_motif3, average_real1_motif4 = find_average(G1, 1)

average_random2_motif3, average_random2_motif4 = find_average(randomG2, 100)
average_real2_motif3, average_real2_motif4 = find_average(G2, 1)

columns_3 = ["motif1", "motif2"]
columns_4 = ["motif1", "motif2", "motif3", "motif4", "motif5", "motif6"]

# display table
df3 = pd.DataFrame([average_real1_motif3, average_random1_motif3, average_real2_motif3, average_random2_motif3],
                   index= ["real 1", "random 1", "real 2", "random 2"],
                   columns= columns_3)

df4 = pd.DataFrame([average_real1_motif4, average_random1_motif4, average_real2_motif4, average_random2_motif4],
                   index= ["real 1", "random 1", "real 2", "random 2"],
                   columns= columns_4)


def save_dataframe_as_table(inputdf, filename):

    ax = plt.subplot(111, frame_on=False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    table(ax, inputdf, loc='center')
    plt.savefig(filename)
    plt.show()


save_dataframe_as_table(df3, "Table_3_motifs.png")
save_dataframe_as_table(df4, "Table_4_motifs.png")


def calculate_z_score(realGraph, randomGraph):
    zscore = (randomGraph - np.mean(realGraph)) / np.std(realGraph)
    return zscore


z_score_1_motif3 = calculate_z_score(average_real1_motif3, average_random1_motif3)
z_score_1_motif4 = calculate_z_score(average_real1_motif4, average_random1_motif4)
z_score_2_motif3 = calculate_z_score(average_real2_motif3, average_random2_motif3)
z_score_2_motif4 = calculate_z_score(average_real2_motif4, average_random2_motif4)

plt.plot(columns_3, z_score_1_motif3)
plt.title("3 node motifs for Random 1 graphs")
plt.savefig('3MotifsRandom1.png')
plt.show()

plt.plot(columns_4, z_score_1_motif4)
plt.title("4 node motifs for Random 1 graphs")
plt.savefig('4MotifsRandom1.png')
plt.show()


plt.plot(columns_3, z_score_2_motif3)
plt.title("3 node motifs for Random 2 graphs")
plt.savefig('3MotifsRandom2.png')
plt.show()


plt.plot(columns_4, z_score_2_motif4)
plt.title("4 node motifs for Random 2 graphs")
plt.savefig('4MotifsRandom2.png')
plt.show()

