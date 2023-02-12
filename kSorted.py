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
        self.p0 = nei0 / (N - 1)
        self.np_random = np.random.default_rng(self.seed)
        self.graph0 = nx.fast_gnp_random_graph(N, self.p0, seed=seed)
        self.kappa_values = []
        self.zone0_nodes = []
        self.zone1_nodes = []
        self.zoneAssigning()
        

    def zoneAssigning(self):
        ktemp = np.array(self.graph0.degree())
            
        kList = ktemp[:,1]
        nodes = ktemp[:,0]
        kListSort = np.sort(kList)
        kLimInd = int(len(kList)/2)
        kLim = kListSort[kLimInd]
        #totalNodes = np.array(graphInitial.nodes())
    
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
        
    def Pk(kList): #k and Pk values 
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
        size_max = np.zeros((len(r0_list), self.N), dtype=np.float64)
        zone_p = np.ones(2)
        for ind_r0, r0 in enumerate(r0_list):
            num_rem = 0
            graph = self.graph0.copy()
            kappa = 100
            zone_p[0] = zone_p[1] * (1 - r0)
            kappaValues = []
            self.kappa_values = []
            while kappa > 2:
                random_node = self.np_random.choice(np.array(graph.nodes()))
                if random_node in self.zone0_nodes:
                    zone_num = 0
                elif random_node in self.zone1_nodes:
                    zone_num = 1
                else:
                    raise ValueError("error")
                p_removal = zone_p[zone_num]
                if self.np_random.random() < p_removal:
                    num_rem += 1
                else:
                    continue
                graph.remove_node(random_node)
                k_temp = list(graph.degree())
                k_list = np.array([i[1] for i in k_temp])
                kM = np.mean(k_list)
                k2List = k_list**2
                k2M = np.mean(k2List)
                kappa = k2M/kM
                kappaValues.append(kappa)
                clustSizes = [len(c) for c in nx.connected_components(graph)]
                relSize = max(clustSizes)/len(graph.nodes())         
                size_max[ind_r0, num_rem-1] = relSize # relative size of the largest cluster 
                
                
                
                