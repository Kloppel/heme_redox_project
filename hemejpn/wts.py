import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, re, subprocess, sys, json, glob
from molmod import *
from molmod.io import FCHKFile
from molmod.io.xyz import XYZReader, XYZFile
import requests

print("""
|__database
|    |__logfiles
|    |  |__*
|    |
|    |__logfilessplit
|    |__pdb
|       |__from_rcsb_old
|       |__prepared
|
|__heme.ipynb
|__split_traj_com.tcl
|__plots
|__tables
""")


erep_pattern = re.compile("nuclear repulsion energy")
dipole_pattern = "Dipole moment (field-independent basis, Debye)"
edisp_pattern = re.compile("Nuclear repulsion after empirical dispersion term")
homo_pattern = re.compile("Alpha  occ. eigenvalues")
polarizability_ex_pattern = re.compile("Dipole polarizability, Alpha")


def TCL_Skript():
    os.chdir("database/pdb")
    os.system('vmd -dispdev text -e split_traj_com.tcl ')
    os.system('mv *.pdb from_rcsb_old')
    os.chdir("../..")


class porphyr:
    defaultdihedlist = [[0, 4, 2, 1], [0, 4, 2, 3], [0, 4, 6, 5], [0, 4, 6, 19], [0, 10, 8, 7], [0, 10, 8, 9], [0, 10, 12, 1], [0, 10, 12, 11], [0, 16, 14, 13], [0, 16, 14, 15], [0, 16, 18, 7], [0, 16, 18, 17], [0, 22, 20, 19], [0, 22, 20, 21], [0, 22, 24, 13], [0, 22, 24, 23], [1, 2, 3, 5], [1, 2, 4, 0], [1, 2, 4, 6], [1, 12, 10, 0], [1, 12, 10, 8], [1, 12, 11, 9], [2, 1, 12, 10], [2, 1, 12, 11], [2, 4, 0, 10], [2, 4, 0, 16], [2, 4, 0, 22], [2, 4, 6, 19], [3, 2, 1, 12], [3, 2, 4, 0], [3, 5, 6, 19], [4, 0, 10, 8], [4, 0, 10, 12], [4, 0, 16, 14], [4, 0, 16, 18], [4, 0, 22, 20], [4, 0, 22, 24], [4, 2, 1, 12], [4, 6, 19, 20], [5, 3, 2, 1], [5, 6, 4, 0], [5, 6, 19, 20], [6, 4, 0, 10], [6, 4, 0, 16], [6, 4, 0, 22], [6, 4, 2, 1], [6, 19, 20, 21], [6, 19, 20, 22], [7, 8, 9, 11], [7, 8, 10, 0], [7, 8, 10, 12], [7, 18, 16, 0], [7, 18, 16, 14], [7, 18, 17, 15], [8, 7, 18, 16], [8, 7, 18, 17], [8, 10, 0, 4], [8, 10, 0, 16], [8, 10, 0, 22], [8, 10, 12, 1], [9, 8, 7, 18], [9, 8, 10, 0], [9, 11, 12, 1], [10, 0, 4, 2], [10, 0, 4, 6], [10, 0, 16, 14], [10, 0, 16, 18], [10, 0, 22, 20], [10, 0, 22, 24], [10, 8, 7, 18], [10, 12, 1, 2], [11, 9, 8, 7], [11, 12, 1, 2], [11, 12, 10, 0], [12, 1, 2, 3], [12, 1, 2, 4], [12, 10, 0, 4], [12, 10, 0, 16], [12, 10, 0, 22], [12, 10, 8, 7], [13, 14, 15, 17], [13, 14, 16, 0], [13, 14, 16, 18], [13, 24, 22, 0], [13, 24, 22, 20], [13, 24, 23, 21], [14, 13, 24, 22], [14, 13, 24, 23], [14, 16, 0, 4], [14, 16, 0, 10], [14, 16, 0, 22], [14, 16, 18, 7], [15, 14, 13, 24], [15, 14, 16, 0], [15, 17, 18, 7], [16, 0, 4, 2], [16, 0, 4, 6], [16, 0, 10, 8], [16, 0, 10, 12], [16, 0, 22, 20], [16, 0, 22, 24], [16, 14, 13, 24], [16, 18, 7, 8], [17, 15, 14, 13], [17, 18, 7, 8], [17, 18, 16, 0], [18, 7, 8, 9], [18, 7, 8, 10], [18, 16, 0, 4], [18, 16, 0, 10], [18, 16, 0, 22], [18, 16, 14, 13], [19, 6, 4, 0], [19, 6, 4, 2], [19, 6, 5, 3], [19, 20, 21, 23], [19, 20, 22, 0], [19, 20, 22, 24], [20, 19, 6, 4], [20, 19, 6, 5], [20, 22, 0, 4], [20, 22, 0, 10], [20, 22, 0, 16], [20, 22, 24, 13], [21, 20, 19, 6], [21, 20, 22, 0], [21, 23, 24, 13], [22, 0, 4, 2], [22, 0, 4, 6], [22, 0, 10, 8], [22, 0, 10, 12], [22, 0, 16, 14], [22, 0, 16, 18], [22, 20, 19, 6], [22, 24, 13, 14], [23, 21, 20, 19], [23, 24, 13, 14], [23, 24, 22, 0], [24, 13, 14, 15], [24, 13, 14, 16], [24, 22, 0, 4], [24, 22, 0, 10], [24, 22, 0, 16], [24, 22, 20, 19]]
    calldict = {}
    def namema(self,liste):
        return [self.calldict[i] for i in liste]

    def getorders(self,lis):
        return [self.order.index(i) for i in lis]

    def getorder(self,elem):
        return self.order.index(elem)

    def sortorder(self,lis):
        return sorted(lis, key = self.getorder)

    def take_first(self,elem):
        return self.order.index(elem[0])

    def take_second(self,elem):
        return self.order.index(elem[1])

    def take_third(self,elem):
        return self.order.index(elem[2])

    def take_fourth(self,elem):
        return self.order.index(elem[3])


#____________________________________________________________________________________ possible tools for debugging  ____________________________________________________

    def get_compassid(self,listofindizies):
        return [self.compassdict[i] for i in listofindizies]

    def compassorder(self):
        return( [self.get_compassid([self.order[i] for i in dfd ]  ) for dfd in self.defaultdihedlist] )     

    def get_compassname(self,listofindizies):
        strlist = [self.compassdict[i] for i in listofindizies]
        name = strlist[0]
        for n  in strlist[1:] :
            name = name +"_"+n #strlist[i]
        return name
    def compassordername(self):
        return( [self.get_compassname([self.order[i] for i in dfd ]  ) for dfd in self.defaultdihedlist] )     

    def neighbourlist(self):
        for i in self.order:
            print(i," - ",list(self.mol.graph.neighbors[i]))

    def porphycompass(self,Fe, NC2O, C1_NO, C3_NO, N_N, C3_NW, C1_NW ,WC2N, C1_WN, C3_WN, N_W, C3_WS, C1_WS, SC2W, C1_SW, C3_SW, N_S, C3_SO, C1_SO , OC2S, C1_OS, C3_OS, N_O, C3_ON, C1_ON  ):
        print("          {}---{}					".format(C3_NW,C3_NO))
        print("          /      \ ")
        print("        {}       {}					".format(C1_NW,C1_NO))
        print("       /  \     /  \ ")
        print("    {}       {}      {}					".format(WC2N,N_N,NC2O))
        print("     |       |        | 					")
        print("  __{}       |       {}__					".format(C1_WN,C1_ON))
        print("{}    \      |      /     {}					".format(C3_WN,C3_ON))
        print("|      {}----{}----{}      |  					".format(N_W,Fe,N_O))
        print("{}__  /      |      \   __{}					".format(C3_WS,C3_OS))
        print("    {}       |       {}-					".format(C1_WS,C1_OS))
        print("     |       |        | 					")
        print("    {}       {}      {}					".format(SC2W,N_S,OC2S))
        print("       \   /    \   /  ")
        print("        {}       {}					".format(C1_SW,C1_SO))
        print("          \     /")
        print("          {}---{}						".format(C3_SW,C3_SO ))


    def porphycompassnumber(self):
        self.porphycompass(self.Fe, self.NC2O, self.C1_NO, self.C3_NO, self.N_N, self.C3_NW, self.C1_NW ,self.WC2N, self.C1_WN, self.C3_WN, self.N_W, self.C3_WS, self.C1_WS, self.SC2W, self.C1_SW, self.C3_SW, self.N_S, self.C3_SO, self.C1_SO , self.OC2S, self.C1_OS, self.C3_OS, self.N_O, self.C3_ON, self.C1_ON  )
        return "Jo"


    def porphyrcompassname(self):
        self.porphycompass(*self.get_compassid([self.Fe, self.NC2O, self.C1_NO, self.C3_NO, self.N_N, self.C3_NW, self.C1_NW , self.WC2N, self.C1_WN, self.C3_WN, self.N_W, self.C3_WS, self.C1_WS, self.SC2W, self.C1_SW, self.C3_SW, self.N_S, self.C3_SO, self.C1_SO , self.OC2S, self.C1_OS, self.C3_OS, self.N_O, self.C3_ON, self.C1_ON]))

#____________________________________________________________________________________ possible tools for debugging  _____________________________________________________

    def set_porphyr(self):
        angles = []
        for i1 in range(self.mol.size):
            n = list(self.mol.graph.neighbors[i1])
        for i1 in range(self.mol.size):
            n = list(self.mol.graph.neighbors[i1])
            if self.mol.numbers[i1] == 26:
                print(f"Eisen {i1}")
                nofn = []
                NDict = {}
                for N in n:
                    C1 = list(self.mol.graph.neighbors[N])
                    C1Dict ={}
                    for nn in n:
                        try:
                            C1.remove(nn)
                        except:
                            2
                    #print("c1",[i for i in C1 if i != i1])
                    for a in [i for i in C1 if i != i1]:
                        print(a, list(self.mol.graph.neighbors[a]))

                    for c1 in C1:
                        C2 = list(self.mol.graph.neighbors[c1])
                        C1Dict.update({c1:C2})
                        for nn in n:
                            try:
                                C2.remove(nn)
                            except:
                                2
                    NDict.update({N:C1Dict})
                allof = [[[c2 for c2 in c1 ]for c1 in N.values()] for N in NDict.values()]
                c =  []
                for a in allof:
                    for b in a:
                        c = c + b
                self.por_index = []
                for k,N in zip (NDict.keys(),NDict.values()):
                    for c1 in N.values():
                        for c2 in c1:
                            cc=0
                            if c.count(c2) == 2:
                                cc+=1
                            if cc > 0:
                                self.por_index.append(k)
                                for zz in c1:
                                    self.por_index.append(zz)
                                for zz in list(N.keys()):
                                    self.por_index.append( zz )
                nofn = list(set(nofn))
            for index, i0 in enumerate(n):
                for i2 in n[:index]:
                    angles.append((i0, i1, i2))
        self.graph  = self.mol.graph.get_subgraph(self.por_index)

        N = []
        C1 = []
        C2 = []
        C3 = []
        for i in set(self.por_index):
            if self.graph.numbers[i] == 7 :
                N.append(i)
            elif 7 in  self.graph.numbers[ list(self.graph.neighbors[i])] and self.graph.numbers[i]!=26 :
                C1.append(i)
        for i in set(self.por_index):
            if i not in N and i not in C1:
                if all( i in C1 for i in list(self.graph.neighbors[i])):
                    C2.append(i)
                elif self.graph.numbers[i]!=26:
                    C3.append(i)
        self.calldict = {}
        for i in set(self.por_index):
            if i in N:
                self.calldict[i] = "N"
            if i in C1:
                self.calldict[i] = "C1"
            if i in C2:
                self.calldict[i] = "C2"
            if i in C3:
                self.calldict[i] = "C3"
            if self.graph.numbers[i]==26:
                self.calldict[i] = "Fe"
        N1  = N[0]
        C11a,C11b = 1,2
        C11a,C11b = [i for i in  list(self.mol.graph.neighbors[N1]) if i in C1]
        C1C2 = [i for i in list(self.mol.graph.neighbors[C11b]) if i in C2][0]
        C4C1 = [i for i in list(self.mol.graph.neighbors[C11a]) if i in C2][0]
        C31b = [i for i in list(self.mol.graph.neighbors[C11b]) if i in C3][0]
        C31a = [i for i in list(self.mol.graph.neighbors[C11a]) if i in C3][0]
        C12a = [i for i in list(self.mol.graph.neighbors[C1C2]) if i in C1 and i!=C11b]  [0]
        N2 = [i for i in list(self.mol.graph.neighbors[C12a]) if i in N]  [0]
        C32a = [i for i in list(self.mol.graph.neighbors[C12a]) if i in C3 ]  [0]
        C32b = [i for i in list(self.mol.graph.neighbors[C32a]) if i in C3]  [0]
        C12b = [i for i in list(self.mol.graph.neighbors[C32b]) if i in C1]  [0]
        C12a = [i for i in list(self.mol.graph.neighbors[C32a]) if i in C1]  [0]
        C2C3 = [i for i in list(self.mol.graph.neighbors[C12b]) if i in C2][0]
        C13a = [i for i in list(self.mol.graph.neighbors[C2C3]) if i in C1 and i != C12b]  [0]
        C33a = [i for i in list(self.mol.graph.neighbors[C13a]) if i in C3]  [0]
        C33b = [i for i in list(self.mol.graph.neighbors[C33a]) if i in C3]  [0]
        C13b = [i for i in list(self.mol.graph.neighbors[C33b]) if i in C1]  [0]
        C13a = [i for i in list(self.mol.graph.neighbors[C33a]) if i in C1]  [0]
        N3 = [i for i in list(self.mol.graph.neighbors[C13a]) if i in N]  [0]
        C3C4 = [i for i in list(self.mol.graph.neighbors[C13b]) if i in C2][0]
        C14a = [i for i in list(self.mol.graph.neighbors[C3C4]) if i in C1 and i != C13b]  [0]
        C34a = [i for i in list(self.mol.graph.neighbors[C14a]) if i in C3 ]  [0]
        C34b = [i for i in list(self.mol.graph.neighbors[C34a]) if i in C3]  [0]
        C14b = [i for i in list(self.mol.graph.neighbors[C34b]) if i in C1]  [0]
        C14a = [i for i in list(self.mol.graph.neighbors[C34a]) if i in C1]  [0]
        N4 = [i for i in list(self.mol.graph.neighbors[C14a]) if i in N]  [0]
        Fe = [i for i in self.calldict.keys() if self.calldict[i] == "Fe"][0]
        print(f"Fe {Fe}")
        if not all([     all(  i in C1 for i in [C11a,C11b,C12b,C12a]),          all(  i in C2 for i in [C1C2,C4C1]),          all(  i in C3 for i in [C31a,C31b,C32a,C32b]),          all(  i in N for i in [N1,N2])  ]):
            raise ValueError('Fehler bei der Zuweisung.')
        methyldict={}
        for i in C3:
            c = [c for c  in list(self.mol.graph.neighbors[i]) if c not in C3 and c not in C1][0]
            if all([i != 6 for i in self.graph.numbers[[v for v in list(self.mol.graph.neighbors[c]) if v != i]  ]]):
                methyldict[i] = True
            else:
                methyldict[i] = False
            list(self.graph.neighbors[i])
        pairs = [(C31b, C32a), (C32b, C33a), (C33b, C34a), (C34b, C31a)]
        npairs = [N1, N2, N3, N4]
        C1pairs = [(C11b, C12a), (C12b, C13a), (C13b, C14a), (C14b, C11a)]

        self.g1 = [C4C1, C11a, C31a, N1, C31b, C11b]
        self.g2 = [C1C2, C12a, C32a, N2, C32b, C12b]
        self.g3 = [C2C3, C13a, C33a, N3, C33b, C13b]
        self.g4 = [C3C4, C14a, C34a, N4, C34b, C14b]

        for p in pairs:
            if methyldict[p[0]] and methyldict[p[1]]:
                methyl = list(p)
                NS = npairs[pairs.index(p)]
                NN = npairs[pairs.index(p)-2]
            if not methyldict[p[0]] and not methyldict[p[1]]:
                acid = list(p)
                NW = npairs[pairs.index(p)]
                NO  = npairs[pairs.index(p)-2]
 
        for i,g in enumerate( [self.g1,self.g2,self.g3,self.g4]):
            if any(m in g for m in methyl):
                if any(a in g for a in acid):
                    SG = g
                else:
                    OG = g
            else:
                if any(a in g for a in acid):
                    WG = g
                else:
                    NG = g
        WC2N, C1_NW, C3_NW, N_N, C3_NO, C1_NO = NG
        NC2O, C1_ON, C3_ON, N_O, C3_OS, C1_OS = OG
        OC2S, C1_SO, C3_SO, N_S, C3_SW, C1_SW = SG
        SC2W, C1_WS, C3_WS, N_W, C3_WN, C1_WN = WG
        
        if not (C1_OS in list(self.mol.graph.neighbors[OC2S]) and C1_SO in list(self.mol.graph.neighbors[OC2S])):
            NC2O, C1_NO, C3_NO, N_N, C3_NW, C1_NW = NG
            WC2N, C1_WN, C3_WN, N_W, C3_WS, C1_WS = OG
            SC2W, C1_SW, C3_SW, N_S, C3_SO, C1_SO = SG
            OC2S, C1_OS, C3_OS, N_O, C3_ON, C1_ON = WG

            WC2N, C1_NW, C3_NW, N_N, C3_NO, C1_NO = NG
            SC2W, C1_WS, C3_WS, N_W, C3_WN, C1_WN = OG
            OC2S, C1_SO, C3_SO, N_S, C3_SW, C1_SW = SG
            NC2O, C1_ON, C3_ON, N_O, C3_OS, C1_OS = WG

        self.Fe, self.NC2O, self.C1_NO, self.C3_NO, self.N_N, self.C3_NW, self.C1_NW , self.WC2N, self.C1_WN, self.C3_WN, self.N_W, self.C3_WS, self.C1_WS, self.SC2W, self.C1_SW, self.C3_SW, self.N_S, self.C3_SO, self.C1_SO , self.OC2S, self.C1_OS, self.C3_OS, self.N_O, self.C3_ON, self.C1_ON  =Fe, NC2O, C1_NO, C3_NO, N_N, C3_NW, C1_NW ,WC2N, C1_WN, C3_WN, N_W, C3_WS, C1_WS, SC2W, C1_SW, C3_SW, N_S, C3_SO, C1_SO , OC2S, C1_OS, C3_OS, N_O, C3_ON, C1_ON

        self.importantdihed = [Fe, NC2O, OC2S, SC2W, WC2N]
        self.order  = [Fe] + NG + OG + SG + WG
        self.order  = [Fe] + [ NC2O, C1_NO, C3_NO, N_N, C3_NW, C1_NW ] + [WC2N, C1_WN, C3_WN, N_W, C3_WS, C1_WS] + [SC2W, C1_SW, C3_SW, N_S, C3_SO, C1_SO ] + [OC2S, C1_OS, C3_OS, N_O, C3_ON, C1_ON] 
        self.saddling = [[C1_SW, N_S ,N_N,C1_NW],[C1_SO, N_S,N_N, C1_NO],[C1_WS, NW, NO ,C1_OS],[C1_WN, NW, NO ,C1_ON]]
        self.ruffling = [[C3_SW, C1_SW, C1_WS, C3_WS],[C3_SO, C1_SO, C1_OS, C3_OS],[C3_WN, C1_WN, C1_NW, C3_NW],[C3_NO, C1_NO, C1_ON, C3_ON]]
        
        self.compassdict = { Fe:"Fe",   NC2O:"NC2O", C1_NO:"C1_NO",C3_NO:"C3_NO", N_N:"N_N", C3_NW:"C3_NW", C1_NW:"C1_NW",             OC2S:"OC2S", C1_OS:"C1_OS", C3_OS:"C3_OS", N_O:"N_O", C3_ON:"C3_ON", C1_ON:"C1_ON",             SC2W:"SC2W", C1_SW:"C1_SW", C3_SW:"C3_SW", N_S:"N_S", C3_SO:"C3_SO", C1_SO:"C1_SO",            WC2N:"WC2N", C1_WN:"C1_WN", C3_WN:"C3_WN", N_W:"N_W", C3_WS:"C3_WS", C1_WS:"C1_WS" }
        for rs in self.order:
            print(rs, list(self.mol.graph.neighbors[rs]))

    def get_dihed_per_list(self,givenlist):
        dihedral = []
        dihedraldict = {}
        dihed_2 = []
        olist = []        
        for i1, i2, i3, i4 in givenlist:
#            i1, i2, i3, i4 = [self.order[i] for i in [io1, io2, io3, io4]]
            dihedral.append([[i1, i2, i3,i4],self.namema([i1, i2, i3,i4]), dihed_angle(self.mol.coordinates[[i1, i2, i3,i4]])[0]/deg , [self.order.index(i) for i in [i1,i2,i3,i4]]   ])
        return dihedral

    def get_dihed(self):
        dihedlist = [[self.order[i] for i in [io1, io2, io3, io4]]  for io1, io2, io3, io4 in self.defaultdihedlist]
        return self.get_dihed_per_list(dihedlist)

    def get_saddling_compass(self):
        return ["C1_SW_N_S_N_N_C1_NW" ,"C1_SO_N_S_N_N_C1_NO", "C1_WS_NW_NO_C1_OS", "C1_WN_NW_NO_C1_ON"]

    def get_ruffling_compass(self):
        return  ["C3_SW_C1_SW_C1_WS_C3_WS","C3_SO_C1_SO_C1_OS_C3_OS","C3_WN_C1_WN_C1_NW_C3_NW","C3_NO_C1_NO_C1_ON_C3_ON"]


    def get_dihed_ruffling(self):
        return self.get_dihed_per_list(self.ruffling)
    def get_dihed_saddling(self):
        return self.get_dihed_per_list(self.saddling)

    def __init__(self, name):
       # print(name)
        self.t = name
        self.mol = Molecule.from_file(name)
        self.mol.set_default_graph()
        self.set_porphyr()

    def get_angleold(self):
        angles = []
        back = []
        for i1 in self.por_index:
            n = list(self.graph.neighbors[i1])
            for index, i0 in enumerate(n):
                for i2 in n[:index]:
                    angles.append((i0, i1, i2))
        for i0, i1, i2 in angles:
            angle = bend_angle(self.mol.coordinates[[i0, i1, i2]])[0]
            back.append(  [[i0, i1, i2], self.namema([i0 , i1 , i2]), angle/deg       ])
        return back

    def __init__(self, name):
        print(name)
        self.saddling_compass = [["C1_SW", "N_S" ,"N_N","C1_NW"],["C1_SO", "N_S","N_N", "C1_NO"],["C1_WS", "NW", "NO" ,"C1_OS"],["C1_WN", "NW", "NO" ,"C1_ON"]]
        self.ruffling_compass = [["C3_SW", "C1_SW", "C1_WS", "C3_WS"],["C3_SO", "C1_SO", "C1_OS", "C3_OS"],["C3_WN", "C1_WN", "C1_NW", "C3_NW"],["C3_NO", "C1_NO", "C1_ON", "C3_ON"]]

        self.mol = Molecule.from_file2(name)
        self.mol.set_default_graph()
        self.set_porphyr()
        self.porphycompassnumber()
        
        self.porphyrcompassname()

class dihedpdb:
    namelist =   ["Fe", "NC2O", "C1_NO ", "C3_NO ", "N_N ", "C3_NW", "C1_NW ", "OC2S", "C1_OS", "C3_OS", "N_O", "C3_ON", "C1_ON ", "SC2W", "C1_SW", "C3_SW", "N_S", "C3_SO", "C1_SO ", "WC2N", "C1_WN", "C3_WN", "N_W", "C3_WS", "C1_WS "]
    df =pd.DataFrame()
    def orderlist(self,liste):
        str = ""
        for i in liste:
            str = str + self.namelist[i]+" "
        return str



    def __init__(self, **kwargs): #, knowndihed=[]):
      
        nr = []
        if kwargs.get("read_keep"):
            try:
                self.df = pd.read_csv("tables/Dihedral.csv")
                self.df = self.df.rename(columns={"Unnamed: 0": "pdb"})
                self.df = self.df.set_index("pdb")
            except:
                self.df = pd.DataFrame()
                print("except_dihed")                

            try:
                self.df_saddling = pd.read_csv("tables/Saddling.csv")
                self.df_saddling = self.df_saddling.rename(columns={"Unnamed: 0": "pdb"})
                self.df_saddling = self.df_saddling.set_index("pdb")
            except:
                self.df_saddling = pd.DataFrame()
                print("except_saddling")  
            try:
                self.df_ruffling = pd.read_csv("tables/Ruffling.csv")
                self.df_ruffling = self.df_ruffling.rename(columns={"Unnamed: 0": "pdb"})
                self.df_ruffling = self.df_ruffling.set_index("pdb")
            except:
                self.df_ruffling = pd.DataFrame()
                print("except ruffling")                                

        else:

            self.df = pd.DataFrame()
            self.df_ruffling = pd.DataFrame()
            self.df_saddling = pd.DataFrame()
       # print(self.df.to_string())
        print(self.df_ruffling.to_string())
        print(self.df_saddling.to_string())
        knowndihed = list(self.df.index)
        knownsaddling = list(self.df_saddling.index)
        knownruffling = list(self.df_ruffling.index)
        for k in glob.glob("database/pdb/prepared/*.pdb"):
            if not (k[k.find("/prepared/")+10:-4] in knowndihed) or not (k[k.find("/prepared/")+10:-4] in knownsaddling) or not (k[k.find("/prepared/")+10:-4] in knownruffling):
                print(k[k.find("/prepared/")+10:-4]," new")
                try:
                    a = porphyr(k)
                    if not (k[k.find("/prepared/")+10:-4] in knowndihed):
                        try:
                            dihed = a.get_dihed()
                            self.df = self.df.append(pd.DataFrame( [[i[2] for i in dihed]], index = [k[k.find("/prepared/")+10:-4]],columns =a.compassordername())) #)   )
                            nr.append([i[3] for i in dihed]) ##überprügen wahrschl müpll
                        except:
                            print(f"problems with {k} dihed" )
                            

                    if not (k[k.find("/prepared/")+10:-4] in knownsaddling):
                        try:
                            saddling = a.get_dihed_saddling()
                           # print([[i[2] for i in dihed]])
                            #print([k[k.find("/prepared/")+10:-4]])
                            #print(a.get_saddling_compass())
                            #print(dihed,"saddling")
                            #print(pd.DataFrame( [[i[2] for i in dihed]], index = [k[k.find("/prepared/")+10:-4]],columns =a.get_saddling_compass()).to_string())
                            print(pd.DataFrame( [[i[2] for i in saddling]], index = [k[k.find("/prepared/")+10:-4]],columns =a.get_saddling_compass()))
                            self.df_saddling = self.df_saddling.append(pd.DataFrame( [[i[2] for i in saddling]], index = [k[k.find("/prepared/")+10:-4]],columns =a.get_saddling_compass()))
                        except:
                            print(f"problems with {k} saddling" )

                    if not (k[k.find("/prepared/")+10:-4] in knownruffling):
                        try:
                            dihed = a.get_dihed_ruffling()
                         #   print(dihed,"ruffling")
                            self.df_ruffling = self.df_ruffling.append(pd.DataFrame( [[i[2] for i in dihed]], index = [k[k.find("/prepared/")+10:-4]],columns =a.get_ruffling_compass()))

                     #       self.df_ruffling = self.df.append(pd.DataFrame( [[i[2] for i in dihed]], index = [k[k.find("/prepared/")+10:-4]],columns =a.ruffling_compass))                         #)   )
                        except:
                            print(f"problems with {k} ruffling" )
                except:
                    print("problems with ",k)
        print(self.df_saddling.to_string())
        
        self.df.index = self.df.index.str.upper()                                
        self.df.to_csv("tables/Dihedral.csv")
        
        self.df_saddling.index = self.df_saddling.index.str.upper()                                
        self.df_saddling.to_csv("tables/Saddling.csv")
        
        self.df_ruffling.index = self.df_ruffling.index.str.upper()                                
        self.df_ruffling.to_csv("tables/Ruffling.csv")


class splitlog_nbo_chloro():
    def __init__():
        """
        The class splitlog_nbo_chloro provides  with the core function.
        """
        return

    def core():
        #

        finaldirectory = "database/logfilessplit/"
        for k in glob.glob('database/logfiles/*.log'):
            with open(k) as f:
                lines = f.readlines()
            indizesofend = [i for i,line  in enumerate(lines) if "Normal termination of Gaussian 16" in line]
            if (len(indizesofend)==2):
                k=k[k.rindex('/')+1:][:-4]
                nbo = lines[:indizesofend[0]+1]
                chloro = lines[indizesofend[0]+1:]
             #   with open(k[k.rindex('/')+1:]+"_nbo.log", 'w') as fp:
                with open( finaldirectory +  k+"_nbo.log", 'w') as fp:                    
                    for line in nbo:
                        fp.write(line)
                with open( finaldirectory +  k+"_chloro.log", 'w') as fp:                    

               # with open(finaldirectory+k[k.rindex('/')+1:]+"_chloro.log", 'w') as fp:
                    for line in chloro:
                        fp.write(line)
            else:
                k=k[k.rindex('/')+1:][:-4]
                with open( finaldirectory +  k+"_nbo.log", 'w') as fp:                    

#                with open(finaldirectory+k+"_nbo.log", 'w') as fp:
                    for line in lines:
                        fp.write(line)
            print("database/logfiles/"+k)
            os.system('mv {} database/logfiles/old/'.format("database/logfiles/"+k+".log"))
        return
#splitlog_nbo_chloro.core()        




class prepare_gaussian_logs:
    finaldirectory = "database/logfilessplit/"
    
    foulder_of_foulders = "database/logfiles/"
    

    def defaultsplit(self,foulder,suffix):
        with open(foulder) as f:
            lines = f.readlines()
        indizesofend = [i for i,line  in enumerate(lines) if "Normal termination of Gaussian 16" in line]
        name=foulder[foulder.rindex('/')+1:][:-4]
        if (name[-2:]) == "05":
            ox = "05"
            pdb = name[:-2]
        elif (name[-2:]) == "12":
            ox = "12"
            pdb = name[:-2]
        elif (name[-2:]) == "16":
            ox = "16"
            pdb = name[:-2]
        elif (name[-2:]) == "01":
            ox = "01"
            pdb = name[:-2]                        
        else:
            ox = "01"
            pdb = name

        if (len(indizesofend)==2):
            
            nbo = lines[:indizesofend[0]+1]
            chloro = lines[indizesofend[0]+1:]
          #  with open( self.finaldirectory +  name+"_nbo.log", 'w') as fp:   
            with open( self.finaldirectory +  pdb+"_"+ ox +"_nbo" + ".log", 'w') as fp:                                    
                for line in nbo:
                    fp.write(line)
#            with open( self.finaldirectory +  name+"_chloro.log", 'w') as fp:                    
            with open( self.finaldirectory +  pdb+"_"+ ox +"_chloro" + ".log", 'w') as fp:                                    

                for line in chloro:
                    fp.write(line)
        else:
            
            with open( self.finaldirectory +  pdb+"_"+ ox +"_nbo" + ".log", 'w') as fp:                                    

                for line in lines:
                    fp.write(line)
        #print("database/logfiles/"+name)

    def core(self):
        for k in glob.glob(self.foulder_of_foulders+"*"):
            j = k[k.rindex("/")+1:]
            l  = glob.glob(k+"/*")
            Linksuffix = {  
                "Link1.log":"_01_nbo.log",
                "Link2.log":"_01_chloro.log",
                "05Link3.log": "_05_nbo.log",
                "05Link4.log":"_05_chloro.log",
                "12Link5.log": "_12_nbo.log",
                "12Link6.log":"_12_chloro.log",
                "16Link7.log": "_16_nbo.log",
                "16Link8.log":"_16_chloro.log"        
            }
            for L in l:
                L2 = L[L.rindex("/")+1:]
               # print(L2)
                try:
                    sfx = Linksuffix[ [i for i in list(Linksuffix.keys()) if i in L2][0]]
                    os.system( "cp {} {}/{}".format(L,self.finaldirectory,j+sfx)     )               
                except:
                    if L[L.rindex("/")+1:-4] == k[k.rindex("/")+1:]:
                        #print(L)
                        print(L,j+"_01.log")
                        self.defaultsplit(L,j+"_01.log")
                    #    print(L,"  ",j+"01.log")

                    else:
                        self.defaultsplit(L,L2)
                       # print(L,"  ",L2)

    def __init__(self):
        return


class helpers():
    def coordtomatrix(coordinates):
        Matrix=[[int(float(i)) if int(float(i))==float(i) else float(i) for i in j.split() ] for j in coordinates]
        df1 =pd.DataFrame(Matrix)
        df = df1.rename({0 : 'CenterNumber', 1: 'AtomicNumber', 2:'AtomicType', 3: 'X' ,4:'Y', 5:'Z'  },axis=1)
        return  df

    def dist(tupel1, tupel2):
        a = [(tupel1[i]-tupel2[i])**2 for i in range(0,3)]
        b=np.sqrt(a[0]+a[1]+a[2])
        return b

    def einmallist(liste):
        LISTE= liste.copy()
        for i in LISTE:
            liste.remove(i)
            if i in liste:
                return False
        return True

    def dreierkombi(zahl):
        listlist = [[str(a),str(b),str(c)] for a in range(1,zahl+1) for b in range(1,zahl+1) for c in range(1,zahl+1) if helpers.einmallist([a,b,c])]
        liste = list(np.array(listlist).reshape(len(listlist)*3))
        return liste

    def viererkombi(zahl):
        listlist = [[str(a),str(b),str(c),str(d)] for a in range(1,zahl+1) for b in range(1,zahl+1) for c in range(1,zahl+1) for d in range(1,zahl+1) if helpers.einmallist([a,b,c,d])]
        liste = list(np.array(listlist).reshape(len(listlist)*4))
        return liste

    def sechserkombi(zahl):
        listlist = [[str(a),str(b),str(c),str(d),str(e),str(f)] for a in range(1,zahl+1) for b in range(1,zahl+1) for c in range(1,zahl+1) for d in range(1,zahl+1) for e in range(1,zahl+1) for f in range(1,zahl+1) if helpers.einmallist([a,b,c,d,e,f])]
        liste = list(np.array(listlist).reshape(len(listlist)*6))
        return liste

    def get_geom(streams): # extracts the geometry from the compressed stream - input orientation!
        geom = []
        for item in streams[-1][16:]:
            if item == "":
                break
            geom.append([item.split(",")[0],float(item.split(",")[-3]),float(item.split(",")[-2]),float(item.split(",")[-1])])
        return(geom)

    def dicttodelete(dictofvalues, number=0):
        arbeitsdict={}
        arbeitsdict.update({
                            'pdb':dictofvalues["pdb"],#[0],
                            'Ox':dictofvalues["Ox"],#[0],
                            'spin':dictofvalues["spin"],#[0], 
                            'method':dictofvalues["method"],#[0],
                            "e":dictofvalues["e"][number] ,
                            "edisp":dictofvalues["edisp"],
                            "homo":dictofvalues["homo"][0],
                            "lumo":dictofvalues["homo"][1],
                            "chem_pot":dictofvalues["homo"][2],
                            "dipole" : dictofvalues["dipole"]  ,
                            "qpole1" : dictofvalues["qpole"][0] ,
                            "qpole2" : dictofvalues["qpole"][1],
                            "qpole3" : dictofvalues["qpole"][2],
                            "qpole4" : dictofvalues["qpole"][3],
                            "polar-iso":  dictofvalues["polar-iso"],
                            "polar-aniso":  dictofvalues["polar-aniso"]})
        return arbeitsdict
    
    def delete():
        def a(pdb):
            dfwork = dfges.loc[[pdb]].copy()
            dfbela = dfwork.copy()
            pdb_I = [i+"_"+c*"I" for i in  dfwork.index for c in range(1,len(dfexp.loc[pdb])+1)]
            for c in range(1,len(dfexp.loc[pdb])):
                dfwork = dfwork.append(dfwork)
            #dfbela
            dfwork["pdb_I"]= pdb_I
            return dfwork


class Physical_quantity:
    

    def get_planeangle(self):     # Das muss effizienter werden
        streams = self.outstreams
        atoms=["2","5","8","12","13","19"]
        atoms = ["1","2","3","4","5","6","7","8","9","10","11", "12"]
        if streams[-1][-1] == "@":
            geom = helpers.get_geom(streams)    # reading .log files
        else:
            geom = streams  # reading .gjf files
        planeangleout = ""
        error = ""
        if len(atoms)%6 != 0:
            error = str(len(atoms)) + " numbers given. Plane angles require sets of 6 numbers each "
        for atom in atoms:
            if not atom.isdigit():
                error += atom + ": Only numbers accepted as input for plane angle "
            if int(atom) > len(geom):
                error += atom + " is out of range. Maximum valid atom number: " + str(len(geom)+1) + " "
        if error != "": return(None,error+";"+int(len(atoms)/6)*";")

        for i in range(int(len(atoms)/6)):
            a = geom[int(atoms[6*i+0])-1][:4] # Atomcoords
            b = geom[int(atoms[6*i+1])-1][:4]
            c = geom[int(atoms[6*i+2])-1][:4]
            d = geom[int(atoms[6*i+3])-1][:4]
            e = geom[int(atoms[6*i+4])-1][:4]
            f = geom[int(atoms[6*i+5])-1][:4]

            ab = np.array([a[1]-b[1],a[2]-b[2],a[3]-b[3]]) # Vectors
            bc = np.array([b[1]-c[1],b[2]-c[2],b[3]-c[3]])
            de = np.array([d[1]-e[1],d[2]-e[2],d[3]-e[3]])
            ef = np.array([e[1]-f[1],e[2]-f[2],e[3]-f[3]])

            n1 = np.cross(ab,bc) # Normal vectors
            n2 = np.cross(de,ef)

            planeangle = round(np.degrees(np.arccos(np.dot(n1,n2) / (np.linalg.norm(n1)*np.linalg.norm(n2)))),3)
            planeangle = min(abs(planeangle),abs(180-planeangle))
            planeangleout += str(a[0]+atoms[4*i+0]+" " + b[0]+atoms[4*i+1]+" " + c[0]+atoms[4*i+2]+" " + d[0]+atoms[4*i+3]) + ";" + str(planeangle) + ";"
        return(planeangleout,error)

    def set_planeangle(self):
        self.planeangle ,self.error["planeangle"]= self.get_planeangle()

    def get_dihedrals(self): # das hier muss auch effizienter werden
        streams = self.outstreams
        atoms = ["2","4","6","8"]
        atoms = [1,2,3,4,5,6,7,8,9,10,11, 12]
        atoms = ["1","2","3","4","5","6","7","8","9","10","11", "12"]
        if streams[-1][-1] == "@":
            geom = helpers.get_geom(streams)    # reading .log files
        else:
            geom = streams  # reading .gjf files
        dihedralout = ""
        error = ""
        # check input
        if len(atoms)%4 != 0:
            error = str(len(atoms)) + " numbers given. Dihedrals require sets of 4 numbers each "
        for atom in atoms:
            if not atom.isdigit():
                error += atom + ": Only numbers accepted as input for dihedrals "
            if int(atom) > len(geom):
                error += atom + " is out of range. Maximum valid atom number: " + str(len(geom)+1) + " "
        if error != "": return(None,error+";"+int(len(atoms)/4)*";")
        for i in range(int(len(atoms)/4)):
            a = geom[int(atoms[4*i+0])-1][:4] # Atomcoords
            b = geom[int(atoms[4*i+1])-1][:4]
            c = geom[int(atoms[4*i+2])-1][:4]
            d = geom[int(atoms[4*i+3])-1][:4]
            ab = np.array([a[1]-b[1],a[2]-b[2],a[3]-b[3]]) # Vectors
            bc = np.array([b[1]-c[1],b[2]-c[2],b[3]-c[3]])
            cd = np.array([c[1]-d[1],c[2]-d[2],c[3]-d[3]])
            n1 = np.cross(ab,bc) # Normal vectors
            n2 = np.cross(bc,cd) - bc
            dihedral = round(np.degrees(np.arccos(np.dot(n1,n2) / (np.linalg.norm(n1)*np.linalg.norm(n2)))),3)
            dihedralout += str(a[0]+atoms[4*i+0]+" " + b[0]+atoms[4*i+1]+" " + c[0]+atoms[4*i+2]+" " + d[0]+atoms[4*i+3]) + ";" + str(dihedral) + ";"
        return(dihedralout,error)

    def set_dihedrals(self):
        self.dihedral ,self.error["dihedral"]= self.get_dihedrals()

    def get_angles(self): # das hier muss auch effizienter werden
        streams = self.outstreams
        if streams[-1][-1] == "@":
            geom = helpers.get_geom(streams)    # reading .log files
        else:
            geom = streams  # reading .gjf files
        anglesout = ""
        error = ""
        atoms = [1,2,3,4,5,6,7,8,9,10,11, 12]
        atoms = ["1","2","3","4","5","6","7","8","9","10","11", "12"]
        if len(atoms)%3 != 0:
            error = str(len(atoms)) + " numbers given. Angles require sets of 3 numbers each "
        for atom in atoms:
            if not atom.isdigit():
                error += atom + ": Only numbers accepted as input for angles "
            if int(atom) > len(geom):
                error += atom + " is out of range. Maximum valid atom number: " + str(len(geom)+1) + " "
        if error != "": return(None,error+int(len(atoms)/3)*";")
        for i in range(int(len(atoms)/3)):
            a = geom[int(atoms[3*i+0])-1][:4] # Atomcoords
            b = geom[int(atoms[3*i+1])-1][:4]
            c = geom[int(atoms[3*i+2])-1][:4]
            ba = np.array(a[1:]) - np.array(b[1:])
            bc = np.array(c[1:]) - np.array(b[1:])
            cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
            angle = np.arccos(cosine_angle)
            anglesout += str(a[0]+atoms[3*i+0]+" " + b[0]+atoms[3*i+1]+" " + c[0]+atoms[3*i+2]) + ";" +str(round(np.degrees(angle),3)) + ";"
        return(anglesout,error)

    def set_angles(self):
        self.angles ,self.error["angle"]= self.get_angles()

    def get_outstreams(self): # gets the compressed stream information at the end of a Gaussian job
        streams = []
        starts,ends = [],[]
        error = "failed or incomplete job" # default unless "normal termination" is in file
        loglines=self.lines
        for i in range(len(loglines)):
            if "1\\1\\" in loglines[i]:
                starts.append(i)
            if "@" in loglines[i]:
                ends.append(i)
            if "Normal termination" in loglines[i]:
                error = ""
        if len(starts) != len(ends) or len(starts) == 0: #probably redundant
            error = "failed or incomplete job"
            return(streams,error)
        for i in range(len(starts)):
            tmp = ""
            for j in range(starts[i],ends[i]+1,1):
                tmp = tmp + loglines[j][1:-1]
            streams.append(tmp.split("\\"))
        return(streams,error)

    def set_outstreams(self):
        self.outstreams=self.get_outstreams()[0]

    def get_distances(self):
        df=self.df.copy()
        distlist = [helper.dist((df["X"][i],df["Y"][i],df["Z"][i]),((df["X"][j],df["Y"][j],df["Z"][j]))) for i in range(0,df.index.stop)  for j in range(0,df.index.stop) ]
        distmatrix = np.array(distlist).reshape(self.df.index.stop,self.df.index.stop)
        distdict = pd.DataFrame(distmatrix)
        return distdict

    def set_distances(self):
        self.distances = self.get_distances()

    def get_edisp(self):
        filecont=self.lines
        error = "no dispersion correction found in this file"
        disp = 0
        for i in range(len(filecont)-1,1,-1):
            if edisp_pattern.search(filecont[i]):
                e_rep_disp = float(filecont[i].split()[-2]) # hartree
                disp = 1
            if erep_pattern.search(filecont[i]) and disp == 1:
                e_rep_nodisp = float(filecont[i].split()[-2]) # hartree
                e_disp = (e_rep_disp - e_rep_nodisp) * 627.50947
                return(e_disp,"")
        return(None,error)

    def set_edisp(self):
        self.edisp, self.error["edisp"]=self.get_edisp()

    def get_polarizability(self): #iso aniso
        filecont = self.lines
        for i in range(len(filecont)-1,1,-1):
            if polarizability_ex_pattern.search(filecont[i]):
                alpha_iso = float(filecont[i+4].split()[1].replace("D","E"))
                alpha_aniso = float(filecont[i+4].split()[2].replace("D","E"))
                #return(str(alpha_iso) + ";" + str(alpha_aniso) + ";","")
                return ({"iso":alpha_iso, "aniso":alpha_aniso},"")
        error = "no polarizability information found;"
        alpha_aniso, alpha_iso = None,None
        return({"iso":alpha_iso, "aniso":alpha_aniso}   ,error)

    def set_polarizability(self):
        self.polar, self.error["polar"] = self.get_polarizability()

    def get_quadrupole(self):
        for item in self.outstreams[-1]:
            if "Quadrupole" in item:
                q = item[11:].split(",")
                q = [float(i) for i in q]
                q_comps = np.array(([q[0],q[3],q[4]],[q[3],q[1],q[5]],[q[4],q[5],q[2]]))
                q_diag = np.linalg.eig(q_comps)[0]
                q_ampl = np.linalg.norm(q_diag)
                results = (np.max(q_diag),-(np.max(q_diag)+np.min(q_diag)),np.min(q_diag),q_ampl)
                return(results,"")
        return((None,None,None,None),"no quadrupole information found")

    def set_quadrupole(self):
        self.qpol, self.error["qpole"] = self.get_quadrupole()

    def get_homolumo(self): # homo,lumo energies and derived values of last job in file
        loglines=self.lines
        error = ""
        for line in loglines[::-1]:
            if  homo_pattern.search(line):
                homo = float(str.split(line)[-1])
                lumo = float(str.split(loglines[loglines.index(line)+1])[4])
                mu =  (homo+lumo)/2 # chemical potential / negative of molecular electronegativity
                eta = lumo-homo     # hardness/softness
                omega = str(round(mu**2/(2*eta),5))  # electrophilicity index
                return((float(homo),float(lumo) ,float(round(mu,5)) ,float(round(eta,5)) ,float(omega)),error)
                return()
        error = "no orbital information found;;"
        return ([None]*5,error )

    def set_homolumo(self):
        self.homo, self.error["homo"] = self.get_homolumo()

    def get_dipole(self,no=0):   #no muss nocheinmal überprüft werden
        filecont = self.lines
        if no != -1:
            for i in range(len(filecont)-1,0,-1):
                if dipole_pattern in filecont[i]:
                    dipole = str.split(filecont[i+1])[-1]
                    return(float(dipole),"")
        if no == -1:
            for i in range(len(filecont)-1):
                if dipole_pattern in filecont[i]:
                    no += 1
                    dipole = str.split(filecont[i+1])[-1]
                    if no == 2: return(float(dipole),"")
        error = "no dipole information found;"
        return (None, error    )

    def set_dipole(self):
        self.dipole, self.error["dipol"] = self.get_dipole()

    def get_DataFrame(self):
        index1 = self.lines.index(' Number     Number       Type             X           Y           Z\n' )+2
        index2 = index1 + int([i for i in self.lines if "NAtoms=" in i][0].split()[1])
        coordinates=self.lines[index1:index2]
        df = helpers.coordtomatrix(coordinates)
        return df

    def set_DataFrame(self):
        self.df = self.get_DataFrame()

    def get_e_hf(self,no):
        stream = self.outstreams[no]
        e_hf = ""
        for item in stream:
            if "HF=" in item:
                e_hf = float(item[3:] )
        return(e_hf,"")

    def set_e(self):
        l = np.array([self.get_e_hf(i) for i in range(0,len(self.outstreams)) ])
        self.e = [float(i) for i in  list(l[:,0])]
        self.error["e"] = list(l[:,1])

    def get_dictionary(self):
        return self.quantitiydict

    def get_error(self):
        return self.error

    def get_bothdicts(self):
        return self.quantitiydict , self.error

    def get_name_ox_method(self):

        return self.name.split("_")

    def set_name_ox_method(self):
        self.pdb, oxspin, self.method = self.get_name_ox_method()
        self.spin,self.ox = {"01":[1,0], "05":[5,0] ,"12":[2,1] ,"16":[6,1]}[oxspin]
        self.error["pdb"], self.error["ox"], self.error["method"] = ["","",""]

    def __init__(self, pwd , name = "1f1f"):
        filename = name +".log"
        filename = pwd
        self.name = filename[filename.rindex("/")+1:][:-4]
        with open(filename) as f:
            self.lines = f.readlines()
        self.error={}
#_______________________set variables_________________________________________________________________
        try:
            self.set_DataFrame()
        except:
            self.df=pd.DataFrame({None:[None]})

        self.set_outstreams()
        self.set_e()
        self.set_edisp()
        self.set_dipole()
        self.set_homolumo()
        self.set_quadrupole()
        #self.set_distances()
        self.set_polarizability()
        self.set_name_ox_method()
      #  self.set_planeangle()
       # self.set_dihedrals()
        #self.set_angles()

#_______________________set variables_________________________________________________________________
#_______________________Make the dict_________________________________________________________________
        self.quantitiydict = {}
        csvdict={}
        # self.quantitiydict["Coordinates"] = self.df.to_numpy()
       # self.quantitiydict["Distances"] = self.distances.to_numpy()
        self.quantitiydict["pdb"] = self.pdb.upper()
        self.quantitiydict["Ox"] = self.ox
        self.quantitiydict["spin"] = self.spin   
        self.quantitiydict["method"] = self.method   
        self.quantitiydict["e"] = [np.mean(self.e)]
        self.quantitiydict["edisp"]   = self.edisp
        self.quantitiydict["homo"] = self.homo
        self.quantitiydict["dipole"] = self.dipole
        self.quantitiydict["qpole"] = self.qpol
        self.quantitiydict["polar-iso"] = self.polar["iso"]
        self.quantitiydict["polar-aniso"] = self.polar["aniso"]
        self.quantitiydict["method"] = self.method   
 

class onecsvold:
    dictoferrors = {}

    csv = "tables/calculated.csv"

    def get_csv(self):
      #  for k in self.dictoflogs.keys():
       #     print(len(self.dictoflogs.values()))
        for k in self.dictoflogs.keys():
            #print(len(self.dictoflogs.values()))            
            dictofvalues=self.dictoflogs[k]
            #print(len(dictofvalues))
            for e,E in enumerate( dictofvalues["e"]):
                if (len(dictofvalues["e"]) ==1):
                    titel = k
                else:
                    titel = k+"_"+{0:"a",1:"b",2:"c",3:"d",4:"e"}[e]
                number = e
                got_answer = False
                writeboolean = True
               # if (titel in self.df.columns):
               #     print(titel," schon enthalten")
                #    while got_answer==False:
                 ##       got_answer = False
                   #     userinput = input("soll das selbige überschrieben werden? [y] ja, [n] nein \n")
                    #    if userinput in ["y","n"]:
                     #       got_answer = True
                    #writeboolean = {"y":True,"n":False}[userinput]
                if(writeboolean):
                    arbeitsdict= helpers.dicttodelete(dictofvalues,number)
                    print(arbeitsdict.values())
                    data=list(arbeitsdict.values())
                    print("arbeitsd : ",list(arbeitsdict.keys()))
                    print("cols  ",self.df.columns)
                    print(self.df.to_string())
                    self.df[titel]=data
        return

    def make_dict(self):
        for k in glob.glob('database/logfilessplit/*.log'):
            #i = k[0:-4]
            i=k[k.rindex('/')+1:][:-4]
            if k in self.df.columns:
                print(k ,"schon erhalten")
            else:
                try:    
                    self.dictoflogs[i], self.dictoferrors[i] = Physical_quantity(k).get_bothdicts()
                    #print(len(self.dictoflogs[i]))
                except:
                    print( i, "failed")



    def save_csv(self):


       # dfcalc1 = pd.read_csv("tables/calculated.csv")
        dft = self.df.T
        header_row = dft.iloc[0]
        df_calc = pd.DataFrame(dft.values[1:], columns=header_row)
        df_calc = df_calc.drop([0]).set_index("pdb")
        df_calc.to_csv("tables/calculated.csv")

        #self.df.to_csv(self.csv)

    def __init__(self):
        self.dictoflogs = {}
        try:
            self.df = pd.read_csv(self.csv, index_col=[0])
        except:
            print("exept")
            newindex = ['pdb','Ox','spin','method',"e","edisp","homo","lumo","chem_pot","dipole" ,"qpole1" ,"qpole2" ,"qpole3" ,"qpole4" ,"polar-iso","polar-aniso"] 
            self.df =pd.DataFrame({ "example":['pdb','Ox','spin','method',"e","edisp","homo","lumo","chem_pot","dipole" ,"qpole1" ,"qpole2" ,"qpole3" ,"qpole4" ,"polar-iso","polar-aniso"] }, index=newindex)
        self.make_dict()
        self.get_csv()
        self.save_csv()
        return


class onecsv:
    dictoferrors = {}

    csv = "tables/calculated.csv"

    def get_csv(self):
        for k in self.dictoflogs.keys():
            dictofvalues=self.dictoflogs[k]
            for e,E in enumerate( dictofvalues["e"]):
                if (len(dictofvalues["e"]) ==1):
                    titel = k
                else:
                    titel = k+"_"+{0:"a",1:"b",2:"c",3:"d",4:"e"}[e]
                number = e
                got_answer = False

                arbeitsdict= helpers.dicttodelete(dictofvalues,number)
                data=list(arbeitsdict.values())
                self.df[titel]=data
        return

    def make_dict(self):
        for k in glob.glob('database/logfilessplit/*.log'):
            #i = k[0:-4]
            i=k[k.rindex('/')+1:][:-4]
            if k in self.df.columns:
                print(k ,"schon erhalten")
            else:
                try:    
                    self.dictoflogs[i], self.dictoferrors[i] = Physical_quantity(k).get_bothdicts()
                    #print(len(self.dictoflogs[i]))
                except:
                    print( i, "failed")



    def save_csv(self):
        dft = self.df.T
        header_row = dft.iloc[0]
        df_calc = pd.DataFrame(dft.values[1:], columns=header_row)
        try:
            df_calc = df_calc.drop([0])
        except:
            2    
        df_calc.set_index("pdb")
        df_calc.to_csv("tables/calculated.csv")

        #self.df.to_csv(self.csv)

    def __init__(self):
        self.dictoflogs = {}

        newindex = ['pdb','Ox','spin','method',"e","edisp","homo","lumo","chem_pot","dipole" ,"qpole1" ,"qpole2" ,"qpole3" ,"qpole4" ,"polar-iso","polar-aniso"] 
        self.df =pd.DataFrame({ "example":['pdb','Ox','spin','method',"e","edisp","homo","lumo","chem_pot","dipole" ,"qpole1" ,"qpole2" ,"qpole3" ,"qpole4" ,"polar-iso","polar-aniso"] }, index=newindex)
        self.make_dict()
        self.get_csv()
        self.save_csv()
        return


def read_redpot_lit():
    dfex = pd.read_csv("tables/RedoxPotential.csv")
    dfex["pH"] = [ ph[1:]  if "~" in ph else ph for ph in dfex["pH"].astype(str) ]
    [ph  for ph in dfex["pH"].astype(str) if "~" in ph ]
    dfex = dfex[["PDB","pH","Elektrode","Ligand","EM","Delta EM","Variant","Einheit"]]
    dfex = dfex.loc[list(dfex["PDB"].dropna().index)]
    dfex["pdb"] = dfex["PDB"].str.lower()
    dfex = dfex.sort_values(by=['pdb'])
    def split_value_std(string):
        splited =  string.replace(" ","").split("+-")
        if len(splited)==2:
            return splited
        else:
            return [splited[0],0]
    dfex["EM"] = dfex["EM"].astype(str)

    vs = [split_value_std(i) for i in dfex["EM"]]
    dfex["EMv"] = [v[0].replace("?","nan") for v in vs]
    dfex["EM_err"] = [s[1] for s in vs]

    dfex["EMv"] = dfex["EMv"].astype(float)
    dfex["EM_err"] = dfex["EM_err"].astype(float)
    dftest = dfex.copy()
    dfex_1 = dftest.drop_duplicates(subset=["pdb"], keep =False)
    dfex_1["pdb"] = dfex_1["pdb"].str.upper()
    dfex_1 = dfex_1.set_index("pdb")
    dfex_2 = dftest[dftest.duplicated(subset="pdb",keep = False)]
    dfex_2.to_csv("tables/duplicated.csv")
    dfex_1.to_csv("tables/mono.csv")



class Hemetype:

    def get_csv(self):
        return self.df

    def save_csv(self):
        self.df.to_csv("tables/Hemetypes.csv")
    
    def make_csv(self):
        self.df = pd.DataFrame({"pdb" : list(self.types.keys()) ,"type": list(self.types.values())       })
        self.df = self.df.set_index("pdb")
        
        return 
    def perapi(self):
        self.types = {}
        for i in self.listofpdb:
            str1 = ""
            r = requests.get('https://www.rcsb.org/structure/'+i, auth=('user', 'pass'))
            if "HEME C" in str(r.content):
                str1 = str1 + "Heme C"
               # self.types.update({i:"Heme C"})
            if "PROTOPORPHYRIN IX" in str(r.content):
                str1 = str1 + "Heme B"

              #  self.types.update({i:"Heme B"})                
            if "HEME-A" in str(r.content):
             #   self.types.update({i:"Heme A"})
                str1 = str1 + "Heme A"

            
            if str1 == "":
                str1 = "unknown"
            self.types.update({i:str1})            

       
        return

    def __init__(self,listofpdb):
        self.listofpdb = listofpdb
        self.perapi()
        self.make_csv()
#a = Hemetype(["2DGE", "1YNR"])
#a.save_csv()