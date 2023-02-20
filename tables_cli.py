import wts
import glob
import pandas as pd
#a = tables.porphyr("")

#pdb = "database/pdb/prepared/1ccg.pdb"
#a = wts.porphyr(pdb)

#wts.dihedpdb()        
wts.prepare_gaussian_logs().core()
wts.onecsv()

#wts.read_redpot_lit()
#wts.TCL_Skript()
#wts.dihedpdb(read_keep=False)
df = pd.read_csv("tables/Ruffling.csv")
l = list(df[ "Unnamed: 0"])
print(l)
wts.Hemetype(l)

