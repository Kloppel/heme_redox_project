class onecsvold:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os, re, subprocess, sys, json, glob
    import molmod
    from molmod.io import FCHKFile
    from molmod.io.xyz import XYZReader, XYZFile
    import requests

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
        return



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