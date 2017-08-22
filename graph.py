# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 12:58:05 2017

@author: Benjamin
"""


"""
>>> n = Node(5)
>>> p = Node(6)
>>> q = Node(7)
>>> n.add_child(p)
>>> n.add_child(q)
>>> n.children
[<__main__.Node object at 0x02877FF0>, <__main__.Node object at 0x02877F90>]
>>> for c in n.children:
...   print c.data
"""

"""
import networkx as nx
import matplotlib.pyplot as plt


def generate_graph(generations):

    G = nx.Graph()
    cntct_chain = []
    for gen in generations:
        cntct_chain.append(gen[0])
        
        
    for bn_idx, bn in enumerate(cntct_chain):
        #print("bn" + str(bn))
        G.add_node(bn_idx, atom_idx = bn[0].tolist())
    
    
    #e=[('a','b',0.3),('b','c',0.9),('a','c',0.5),('c','d',1.2)]
    #G.add_weighted_edges_from(e)
    return G




def plot_graph(G):
    
    nodes = G.nodes()
    print(nodes)
    pos=nx.spring_layout(G)
    nx.draw_networkx(G, pos=pos, node_size=100, font_size = 14, with_labels = True)
    labels={}
    for node in nodes:
        labels[node] = str(G.node[node]["atom_idx"])
    for p in pos:
        pos[p] = [pos[p][0],pos[p][1]-0.1]
    nx.draw_networkx_labels(G,pos,labels,font_size=16, font_color = "b")
    
    plt.axis('off')
    plt.savefig("labels_and_colors.png") # save as png
    plt.show() # display
"""