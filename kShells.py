# -*- coding: utf-8 -*-
"""
Created on Sat Feb  4 15:07:59 2023

@author: argdi
"""

import networkx as nx
import numpy as np


class PercolationModel:
    def __init__(self, N, nei0, seed):
        self.N = N
        self.nei0 = nei0
        self.p0 = nei0 / (N - 1) #probability if a link creation
        self.np_random = np.random.default_rng(seed) 
        self.graph0 = nx.fast_gnp_random_graph(N, self.p0, seed=seed) #Erdos Renyi network
        self.kappa_values = [] # kappa = <k^2> / <k> where <k> is the average degree
        self.shells = [] 
        self._k_shell_decomposition()
        self.zone0_nodes = self.shells[-1] #internal zone
        self.zone1_nodes = [node for shell in self.shells[:-1] for node in shell] #external zone 

    def _k_shell_decomposition(self):
        k = 1
        graph = self.graph0.copy()
        shells_temp = []
 
        while True:
            while True:
                if not graph:
                    break
                k_temp = np.array(graph.degree(), dtype=np.int32)
                if len(k_temp) < 1:
                    break
                node_list = k_temp[:, 0] #nodes of the network
                k_list = k_temp[:, 1] #degree connectivity of the network
                
                #find and remove all the nodes that have degree <= k and continue until you cannot remove nodes anymore
                rem_nodes_ind = np.where(k_list <= k)[0] 
                rem_nodes = node_list[rem_nodes_ind] 
                shells_temp.extend(rem_nodes)
                if rem_nodes:
                    graph.remove_nodes_from(rem_nodes)
                else:
                    #if there are not other nodes left to remove, increase k
                    self.shells.append(shells_temp)
                    k += 1
                    shells_temp = []
            break
        self.shells.append(shells_temp)
        
    def Pk(kList): #calculate the degree distribution, providing the degree values
        d = {}
        kAll = []
        PkAll = []
        for i in range(len(kList)):    
            if kList[i] not in d:
                d[kList[i]] = 1
            else:        
                d[kList[i]] += 1
        Pktemp = np.array(list(d.values()))
        PkAll.append(Pktemp/sum(Pktemp))
        kAll.append(np.array(list(d.keys())))


        return [PkAll, kAll]

    def simulate(self, r0_list):
        #r0_list is the list of the parameter values, r0, that defines the internal zone removal probability 
        size_max = np.zeros((len(r0_list), self.N), dtype=np.float64) #relative size of the largest cluster during the removal process
        zone_p = np.ones(2) #removal probabilities
        for ind_r0, r0 in enumerate(r0_list):
            num_rem = 0 
            graph = self.graph0.copy()
            kappa = 100 # kappa = <k^2> / <k> where <k> is the average degree
            zone_p[0] = zone_p[1] * (1 - r0) # internal zone removal probability, using the parameter r0
            kappaValues = [] 
            #self.kappa_values = []
            while kappa > 2: # kappa = 2 is the network breakdown point
                while True:
                    random_node = self.np_random.choice(np.array(graph.nodes()))
                    if random_node in self.zone0_nodes: # internal zone
                        zone_num = 0
                    elif random_node in self.zone1_nodes: # external zone
                        zone_num = 1
                    else:
                        raise ValueError("error")
                    p_removal = zone_p[zone_num] 
                    if self.np_random.random() < p_removal: # if a random number 0 - 1 is smaller than the removal probability we remove the node
                        num_rem += 1
                        break
                graph.remove_node(random_node) 
                k_temp = list(graph.degree()) 
                k_list = np.array([i[1] for i in k_temp]) #degrees of the nodes 
                
                #we calculate kappa
                kM = np.mean(k_list)
                k2List = k_list**2
                k2M = np.mean(k2List)
                kappa = k2M/kM
                kappaValues.append(kappa)
                clustSizes = [len(c) for c in nx.connected_components(graph)] #size of the network's clusters
                relSize = max(clustSizes)/len(graph.nodes()) #relative size of the maximum cluster
                size_max[ind_r0, num_rem-1] = relSize  
                
                
                
                