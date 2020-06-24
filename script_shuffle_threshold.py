#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 16:41:28 2019

@author: eris
"""

import networkx as nx
from networkx.algorithms import bipartite
import csv
import random
import collections
import pandas as pd
import sys

def add_and_remove_edges(G, bp):    
    '''    
    for each node,    
      add a new connection to random other node, with prob p_new_connection,    
      remove a connection, with prob p_remove_connection    

    operates on G in-place    
    '''                
    new_edges = []    
    rem_edges = []   
    #aux_degree = G.degree()

    for node in G.nodes():
        #print(node)
        # find the other nodes this one is connected to
        val_node_bp = list(G.node[node].values())[0]
        
        if val_node_bp == bp:
            connected = [to for (fr, to) in G.edges(node)]
            
            #print(connected)
            # and find the remainder of nodes, which are candidates for new edges   
            unconnected = [n for n in G.nodes() if not n in connected and val_node_bp != list(G.node[n].values())[0]]
            
            k_degree = len(connected)
            #print("\tdegree:\t {} -- {}".format(node, k_degree))
            
            if len(connected):
                for i in range(k_degree):
                    remove = connected[i]
                    G.remove_edge(node, remove)
                    #print("\tedge removed:\t {} -- {}".format(node, remove))
                    rem_edges.append( (node, remove) )
                connected = []
                    
                    
            if len(unconnected):
                for i in range(k_degree):
                    new = random.choice(unconnected)
                    G.add_edge(node, new)
                    #print("\tnew edge:\t {} -- {}".format(node, new))
                    new_edges.append( (node, new) )
                    unconnected.remove(new)
                    connected.append(new)
                    
    return rem_edges, new_edges

def threshold(G, bp):
    degX,degY=bipartite.degrees(G,nodes_0)
    degATC = dict(degX).values()
    degCIE = dict(degY).values()
    counterATC = collections.Counter(degATC)
    counterCIE = collections.Counter(degCIE)
    c_list = []
    nc_list = []
    nuc_list = []
    if bp == 0:
        for th in sorted(list(counterCIE.keys())):
            #th = 1
            H = nx.Graph()
            #for v in G.nodes(data = True):
            #    if v[1]['bipartite'] == 0:
            #        H.add_node(v[0])
        
            for n in G.nodes(data=True):
                if n[1]['bipartite'] == 0:
                    sourceNode = n[0]
                    s_neighbors = set(G.neighbors(n[0]))
                    for m in G.nodes(data = True):
                        if m[1]['bipartite'] == 0: #### Change to 1 to change the projection to active ingredient
                            targetNode = m[0]
                            t_neighbors = set(G.neighbors(m[0]))
                            if sourceNode != targetNode:
                                if len(s_neighbors & t_neighbors) >= th:
                                    H.add_node(sourceNode)
                                    H.add_node(targetNode)
                                    H.add_edge(sourceNode,targetNode)                  
            components = sorted(nx.connected_components(H), key=len, reverse=True)
            #sum(list(map(lambda c: len(c), components)))
            c_list.append(len(components))
            nodes_connected = sum(list(map(lambda c: len(c), components)))
            nc_list.append(nodes_connected)
            nuc_list.append(len(nodes_0) - nodes_connected)
            #nx.write_graphml(H,'proCIE_th_'+str(th)+'.graphml')
    else:
        
        for th in sorted(list(counterATC.keys())):
            #th = 136
            H = nx.Graph()
            #for v in G.nodes(data = True):
            #    if v[1]['bipartite'] == 1:
            #        H.add_node(v[0])
        
            for n in G.nodes(data=True):
                if n[1]['bipartite'] == 1:
                    sourceNode = n[0]
                    s_neighbors = set(G.neighbors(n[0]))
                    for m in G.nodes(data = True):
                        if m[1]['bipartite'] == 1: #### Change to 1 to change the projection to active ingredient
                            targetNode = m[0]
                            t_neighbors = set(G.neighbors(m[0]))
                            if sourceNode != targetNode:
                                if len(s_neighbors & t_neighbors) >= th:
                                    #print(len(s_neighbors & t_neighbors))
                                    #print(sourceNode + " " + targetNode)
                                    H.add_node(sourceNode)
                                    H.add_node(targetNode)
                                    H.add_edge(sourceNode,targetNode)
                                    
            components = sorted(nx.connected_components(H), key=len, reverse=True)
            c_list.append(len(components))
            nodes_connected = sum(list(map(lambda c: len(c), components)))
            nc_list.append(nodes_connected)
            nuc_list.append(len(nodes_1) - nodes_connected)
        #nx.write_graphml(H,'proATC_th_'+str(th)+'.graphml')
    #degXH,degYH=bipartite.degrees(H,nodes_0)
    #degATCH = dict(degXH).values()
    #degCIEH = dict(degYH).values()
    #counterATCH = collections.Counter(degATCH)
    #counterCIEH = collections.Counter(degCIEH)
    return c_list, nc_list, nuc_list, counterATC, counterCIE

if __name__ == '__main__':
    bp = int(sys.argv[1])
    print("Lectura de archivo: ")
    vdmdata_reduce = pd.read_csv('vdmdata_reduce.csv')
    
    print("Identificación de nodos ...")
    nodes_0 = []
    nodes_1 = []
    for m in vdmdata_reduce.iterrows():
        nodes_0.append(m[1][0]) #ICD
        nodes_1.append(m[1][1]) #ATC
        
    nodes_0 = list(dict.fromkeys(nodes_0))
    nodes_1 = list(dict.fromkeys(nodes_1))
    
    for i in range(100):
        # Build a bipartite graph:
        print("Iteración: "+str(i))
        print("Construcción de grafo ...")
        G = nx.Graph()
        G.add_nodes_from(nodes_0, bipartite=0) # Add the node attribute “bipartite” disease
        G.add_nodes_from(nodes_1, bipartite=1) # active substance
        
        for m in vdmdata_reduce.iterrows():
            enfermedad = m[1][0];
            sustancia = m[1][1];
            G.add_edge(enfermedad, sustancia)
        
        print("Shuffle enlaces ...")
        rem_edges, new_edges = add_and_remove_edges(G, bp)
        print("Threshold del grafo shuffle ...")
        cc, nc, nuc, ca, cd = threshold(G, bp)
        cd = sorted(list(cd.keys()))
        ca = sorted(list(ca.keys()))
        print("Creación de dataframes ...")
        if i == 0:
            dcc = pd.DataFrame(cc)
            dnc = pd.DataFrame(nc)
            dnuc = pd.DataFrame(nuc)
            if bp == 0:
                dcd = pd.DataFrame(cd)
            else:
                dca = pd.DataFrame(ca)
        else: 
            dcc.insert(i, i, cc, True)
            dnc.insert(i, i, nc, True)
            dnuc.insert(i, i, nuc, True)
            if bp == 0:
                dcd.insert(i, i, cd, True)
            else:
                dca.insert(i, i, ca, True)
        
        print(str(len(cc)) + " " + str(len(nc)) + " " + str(len(nuc)) + " " + str(len(ca)) + " " + str(len(cd)))
    
    print("Exportación a archivos CSV de los dataframes ...")
    dcc.to_csv (r'export_dcc_'+str(bp)+'.csv', index = None, header=True)
    dnc.to_csv (r'export_dnc_'+str(bp)+'.csv', index = None, header=True)
    dnuc.to_csv (r'export_dnuc_'+str(bp)+'.csv', index = None, header=True)
    if bp == 0:
        dcd.to_csv (r'export_dcd_'+str(bp)+'.csv', index = None, header=True)
    else:
        dca.to_csv (r'export_dca_'+str(bp)+'.csv', index = None, header=True)
    
        
