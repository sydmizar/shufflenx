# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 21:10:18 2020

@author: BALAMLAPTOP2
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 18:21:07 2020

@author: BALAMLAPTOP2
Input variables: 
    type_proj = ['icd', 'atc']
"""

import networkx as nx
from networkx.algorithms import bipartite
import collections
import pandas as pd
import multiprocessing
import sys
import time
import random
import copy 

#def add_and_remove_edges(C, type_proj, degree_atc, degree_icd):
#    C_shuffled = copy.copy(C)
#    
#    for node in C_shuffled.nodes():
#        print(node)
#        val_node_bp = list(C_shuffled.node[node].values())[0]
#        
#        if val_node_bp == 0:
#            connected = [to for (fr, to) in C_shuffled.edges(node)] 
#            unconnected = [n for n in C_shuffled.nodes() if not n in connected and val_node_bp != list(C_shuffled.node[n].values())[0]]
#            k_degree = len(connected)
#
#            full = connected + unconnected
#            full = sorted(full)
#            
#            if len(connected):
#                for i in range(k_degree):
#                    remove = connected[i]
#                    C_shuffled.remove_edge(node, remove)
#                connected = []
#                    
#                    
#            if len(unconnected):
#                for i in range(k_degree):
#                    new = random.choice(full)
#                    ### Recurrir a valor auxiliar para generar loop dependiendo de los compentarios del doctor.
#                    while degY[new] != 0:
#                        if degY[new] != 0:
#                            degY[new] = degY[new] - 1
# 
#                    C_shuffled.add_edge(node, new)
#                    full.remove(new)
#                    connected.append(new)
#                    
#    return C_shuffled
def gf_evaluate( g, x ): 
    """Evaluate the generating function. g: the generating function as a sequence x: the parameter returns: the value"""
    v = 0
    xx = 1
    for k in range(len(g)):
        v = v + g[k] * xx
        xx = xx * x
    return v
 
def gf_differentiate( g ):
    """Differentiate the generating function. g: the generating function returns: the derived function"""
    dg = g[:]
    for k in range(len(g)):
        dg[k] = g[k] * k
    return dg[1:]

def gf_multiplyx( g ):
    """Multiply the generating function through by x. g: the generating function returns: the new function"""
    return [ 0 ] + g

def suffle_edges(C, degree_atc, degree_icd):
    print("Shuffling edges ... ")
    C_shuffled = copy.copy(C)
    k=0
    iter=2*C_shuffled.size()
    while k<iter:
        r1=random.choice(degree_icd)
        d1=random.choice(list(C_shuffled.neighbors(r1)))
        
        r2=random.choice(degree_icd)
        d2=random.choice(list(C_shuffled.neighbors(r2)))
        
        if (C_shuffled.has_edge(r1,d2)==False) & (C_shuffled.has_edge(r2,d1)==False):
            C_shuffled.add_edge(r1,d2)
            C_shuffled.remove_edge(r1,d1)       
            C_shuffled.add_edge(r2,d1)
            C_shuffled.remove_edge(r2,d2)
            k=k+1
    return C_shuffled

def threshold_analysis(C, th, type_proj, nn, iteration):
    print("Creating new graph by given threshold ... "+str(th))
    H = nx.Graph()
    index_source = 1
    index_target = 1
    for n in C.nodes(data=True):
        if n[1]['bipartite'] == type_proj:
            #print(index_source)
            sourceNode = n[0]
            s_neighbors = set(C.neighbors(n[0]))
            for m in C.nodes(data = True):
                if m[1]['bipartite'] == type_proj: #### Change to 1 to change the projection to active ingredient
                    targetNode = m[0]
                    t_neighbors = set(C.neighbors(m[0]))
                    if sourceNode != targetNode and index_target > index_source:
                        if len(s_neighbors & t_neighbors) >= th:
                            H.add_node(sourceNode)
                            H.add_node(targetNode)
                            H.add_edge(sourceNode,targetNode)
                    index_target += 1
            index_target = 1
            index_source += 1        					
    components = sorted(nx.connected_components(H), key=len, reverse=True)
    nodes_connected = sum(list(map(lambda c: len(c), components)))
    nodes_unconnected = nn - nodes_connected
    lcs = len(components[0])
    degrees = H.degree()
    sum_of_edges = sum(list(dict(degrees).values()))
    avg_degree = sum_of_edges / H.number_of_nodes()
    
    if th == 1:
        degrees = H.degree()
        degree_values = sorted(set(dict(degrees).values()))
        pk = [list(dict(degrees).values()).count(i)/float(nx.number_of_nodes(H)) for i in degree_values]
#        print ("mean degree of g0 = {m}".format(m = gf_evaluate(gf_multiplyx(gf_differentiate(pk)), 1)))
#        print ("variance of g0 = {var}".format(var = gf_evaluate(gf_multiplyx(gf_differentiate(pk)), 2)))
        m = gf_evaluate(gf_multiplyx(gf_differentiate(pk)), 1)
        var = gf_evaluate(gf_multiplyx(gf_differentiate(pk)), 2)
        pc = m/(var - m)
        df_degrees = pd.DataFrame(dict(degrees).items(), columns=['code', 'degree'])
        df_degrees.to_csv('degrees_'+str(th)+'_'+str(type_proj)+'.csv', index = False, encoding = 'utf-8')
        with open("mean_variance_data.txt", "a+") as f:
            f.write(str(iteration)+","+str(th)+","+str(m)+","+str(var)+","+str(pc)+"\n")


#        pd.DataFrame(dict(degX_sh).items(), columns=['code', 'degree'])
#        degX_sh,degY_sh=bipartite.degrees(C_shuffled,nodes_0_c)
#        degATC_sh = dict(degX_sh).values()
#        degCIE_sh = dict(degY_sh).values()
#        counterATC_sh = collections.Counter(degATC_sh)
#        counterCIE_sh = collections.Counter(degCIE_sh)
    
    print("Saving values for the given threshold ..."+str(th))
    if type_proj == 0:
        with open("threshold_shuffle_icd_"+str(iteration)+".txt", "a+") as f:
            f.write(str(iteration)+","+str(th)+","+str(len(components))+","+str(nodes_connected)+","+str(nodes_unconnected)+","+str(lcs)+","+str(avg_degree)+"\n")
#        nx.write_graphml(H,'ICD/projICD_th_'+str(th)+'_'+str(iteration)+'.graphml')
    elif type_proj == 1:
        with open("threshold_shuffle_atc_"+str(iteration)+".txt", "a+") as f:
            f.write(str(iteration)+","+str(th)+","+str(len(components))+","+str(nodes_connected)+","+str(nodes_unconnected)+","+str(lcs)+","+str(avg_degree)+"\n")
#        nx.write_graphml(H,'ATC/projATC_th_'+str(th)+'_'+str(iteration)+'.graphml')
    else:
        print("The option doesn't exist. Try again.")
    

if __name__ == '__main__':

#    type_proj = int(sys.argv[1])
    print("Reading file ...")

    vdmdata = pd.read_csv('vdmdata_reduce.csv', encoding = 'utf-8-sig')
    
    nodes_0 = []
    nodes_1 = []
    for m in vdmdata.iterrows():
        nodes_0.append(m[1][0]) #ICD
        nodes_1.append(m[1][1]) #ATC
        
    nodes_0 = list(dict.fromkeys(nodes_0))
    nodes_1 = list(dict.fromkeys(nodes_1))
    print("Building a bipartite graph ...")
    # Build a bipartite graph:
    G = nx.Graph()
    # Add nodes ATC - ICD
    G.add_nodes_from(nodes_0, bipartite=0) # Add the node attribute “bipartite” disease
    G.add_nodes_from(nodes_1, bipartite=1) # active substance
    
    # Add edges without weight
    for m in vdmdata.iterrows():
        enfermedad = m[1][0];
        #peso = m[1][3];
        sustancia = m[1][1];
        G.add_edge(enfermedad, sustancia)
    
    # Get the largest component
    print("Getting largest component ...")
    components = sorted(nx.connected_components(G), key=len, reverse=True)
    largest_component = components[0]
    C = G.subgraph(largest_component)
        
    degX,degY=bipartite.degrees(C,nodes_0)
    degATC = dict(degX).values()
    degCIE = dict(degY).values()
    counterATC = collections.Counter(degATC)
    counterCIE = collections.Counter(degCIE)
#    nx.degree_histogram(G)
#    GPCIE = bipartite.projected_graph(G, nodes_0)
#    GPATC = bipartite.projected_graph(G, nodes_1)
#    degree_freq = nx.degree_histogram(G)
#    degrees = range(len(degree_freq))
#    plt.figure(figsize=(12, 8)) 
#    #plt.loglog(degrees[m:], degree_freq[m:],'go-') 
#    plt.xlabel('Degree')
#    plt.ylabel('Frequency')
    nodes_0_c = []
    nodes_1_c = []
    for n in C.nodes(data=True):
        if n[1]['bipartite'] == 0:
            nodes_0_c.append(n[0])
        if n[1]['bipartite'] == 1:
            nodes_1_c.append(n[0])
#            
#    A = bipartite.biadjacency_matrix(C, nodes_0_c, nodes_1_c)
    
    #use one less process to be a little more stable
    #p = multiprocessing.Pool(processes = multiprocessing.cpu_count()-5)
    p = multiprocessing.Pool()
    #timing it...
    start = time.time()
#    for file in names:
#        p.apply_async(multip, [file])

    for i in range(50):
        print("Shuffle edges ... iteration "+str(i))
        #H = add_and_remove_edges(C, type_proj, dict(degX), dict(degY))
        C_shuffled = suffle_edges(C, sorted(dict(degX).keys()), sorted(dict(degY).keys()))
        degX_sh,degY_sh=bipartite.degrees(C_shuffled,nodes_0_c)
        degATC_sh = dict(degX_sh).values()
        degCIE_sh = dict(degY_sh).values()
        counterATC_sh = collections.Counter(degATC_sh)
        counterCIE_sh = collections.Counter(degCIE_sh)
        nx.write_graphml(C_shuffled,'networks/bipartite_sh_'+str(i)+'.graphml')
        
        print("Apply threshold analysis to shuffled graph ... "+str(i))
        for th_icd in sorted(list(counterCIE_sh.keys())):
            p.apply_async(threshold_analysis, [C_shuffled, th_icd, 0, len(degY_sh), i])
            
        for th_atc in sorted(list(counterATC_sh.keys())):
            p.apply_async(threshold_analysis, [C_shuffled, th_atc, 1, len(degX_sh), i])
            
#        if type_proj == 0:
#            for th in sorted(list(counterCIE_sh.keys())):
#                p.apply_async(threshold_analysis, [C_shuffled, th, type_proj, len(degY_sh), i])
#        elif type_proj == 1:
#            for th in sorted(list(counterATC_sh.keys())):
#                p.apply_async(threshold_analysis, [C_shuffled, th, type_proj, len(degX_sh), i])

    p.close()
    p.join()
    print("Complete")
    end = time.time()
    print('total time (s)= ' + str(end-start))

