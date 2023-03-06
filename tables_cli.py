import wts
import glob
import pandas as pd
#a = tables.porphyr("")

#pdb = "database/pdb/prepared/1ccg.pdb"
#a = wts.porphyr(pdb)

#wts.dihedpdb()        
#wts.prepare_gaussian_logs().core()
#wts.onecsv()

#wts.read_redpot_lit()
#wts.TCL_Skript()
#wts.dihedpdb(read_keep=False)
#df = pd.read_csv("tables/Ruffling.csv")
#l = list(df[ "Unnamed: 0"])

l = ['3CCX',  '1EWH',  '1YCC',  '6U97',  '1CCE',  '1IQC',  '1DRM',  '1F1F',  '4K8F',  '1IQC',  '1H9X',  '2OYY',  '451C',  '1BE3',  '1MZ4',  '1CO6',  '1A6K',  '1CYO',  '1AW3',  '1LS9',  '2PCB',  '2V08',  '1EHB',  '351C',  '1JP9',  '2C1D',  '1B5A',  '5MCS',  '1CCC',  '1C2R',  '1A6M',  '1CXC',  '2PCB',  '1CCK',  '1CTJ',  '1JBQ',  '1HRC',  '1QPU',  '1MOH',  '1GU2',  '1F1C',  '1JP6',  '1YNR',  '1CRC',  '2GSM',  '1CNO',  '1GIW',  '4MU8',  '2C1D',  '1FJ0',  '1A6G',  '2HPD',  '1EUE',  '1JDL',  '1CMT',  '1EBT',  '1HJT',  '1CCG',  '2GSM',  '1VB6',  '1CRY',  '1DJ1',  '3CU4',  '1FLP',  '1CYI',  '1IO3',  '1ES1',  '1CYJ',  '1C75',  '1FHF',  '1CMQ']
#print(l)
wts.Hemetype(l)

#for k in glob.glob("database/pdb/prepared/*.pdb"):
#    print()