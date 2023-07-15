# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 17:20:33 2023

@author: argdi
"""


import networkx as nx
import numpy as np


class PercolationModel:
    def __init__(self, N, nei0, seed):
        self.N = N
        self.seed = seed
        self.nei0 = nei0
        self.p0 = nei0 / (N - 1) #probability if a link creation
        self.np_random = np.random.default_rng(self.seed)
        self.graph0 = nx.fast_gnp_random_graph(N, self.p0, seed=seed) #Erdos Renyi random network
        self.kappa_values = [] # kappa = <k^2> / <k> where <k> is the average degree
        self.zone0_nodes = [] #internal zone
        self.zone1_nodes = [] #external zone
        self.zoneAssigning()
        

    def zoneAssigning(self):
        ktemp = np.array(self.graph0.degree())
        
        #we sort the nodes based on their degree connectivity
        kList = ktemp[:,1]
        nodes = ktemp[:,0]
        kListSort = np.sort(kList)
        kLimInd = int(len(kList)/2)
        kLim = kListSort[kLimInd]
        
        #then we assign half of the highest degree nodes to the internal zone and the rest to the external
        if (self.seed < 501): # each zone will contain the kLim once for half of the realizations
            zone1Ind = np.where(kList<=kLim)[0]
            self.zone1_nodes = nodes[zone1Ind]
            zone0Ind = np.where(kList>kLim)[0]
            self.zone0_nodes = nodes[zone0Ind]
        else:
            zone1Ind = np.where(kList<kLim)[0]
            self.zone1_nodes = nodes[zone1Ind]
            zone0Ind = np.where(kList>=kLim)[0]
            self.zone0_nodes = nodes[zone0Ind]
        
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
            kappa = 100 #kappa = <k^2> / <k> where <k> is the average degree
            zone_p[0] = zone_p[1] * (1 - r0) #internal zone removal probability
            kappaValues = []
            #self.kappa_values = []
            while kappa > 2: # kappa = 2 is the network breakdown point
                while True:
                    random_node = self.np_random.choice(np.array(graph.nodes()))
                    if random_node in self.zone0_nodes: #internal zone
                        zone_num = 0
                    elif random_node in self.zone1_nodes: #external zone
                        zone_num = 1
                    else:
                        raise ValueError("error")
                    p_removal = zone_p[zone_num]
                    if self.np_random.random() < p_removal: #if a random number 0 - 1 is smaller than the removal probability we remove the node
                        num_rem += 1
                        break
                graph.remove_node(random_node)
                k_temp = list(graph.degree())
                k_list = np.array([i[1] for i in k_temp])  #degrees of the nodes 
                
                #we calculate kappa
                kM = np.mean(k_list)
                k2List = k_list**2
                k2M = np.mean(k2List)
                kappa = k2M/kM
                kappaValues.append(kappa)
                clustSizes = [len(c) for c in nx.connected_components(graph)] #size of the network's clusters
                relSize = max(clustSizes)/len(graph.nodes()) #relative size of the largest cluster 
                size_max[ind_r0, num_rem-1] = relSize 
                
                
                
                