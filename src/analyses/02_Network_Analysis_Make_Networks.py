import numpy as np
import pandas as pd
import glob
import os
import pickle
import time
import json

def load_convert_dict():
    with open(fp_net+'Homo_sapiens__ENSP-to-ENSG__All-Mappings.pickle', 'rb') as mydict:
        convert_dict = pickle.load(mydict)
    return convert_dict

def main(File_org):
    Com_dict = {}
    Exp_dict = {}
    with open(File_org, mode='r') as f_old:
        for idx, line in enumerate(f_old):
            if idx == 0:
                continue
            else:
                # get the uniprot IDs and final and initial scores
                ENSPa = line.split(' ')[0].split('.')[1]
                ENSPb = line.split(' ')[1].split('.')[1]
                ExpW = int(line.split(' ')[6])
                ComW = int(line.split(' ')[9])
            # convert ENSP to entrez
            try:
                EntAs = convert_dict[ENSPa]
                EntBs = convert_dict[ENSPb]
                for EntA in EntAs:
                    for EntB in EntBs:
                        EntA = EntA
                        EntB = EntB
                        # sort EntA and EntB and remove self edges
                        if EntA == EntB:
                            continue
                        elif EntB < EntA:
                            EntB, EntA = EntA, EntB
                        else:
                            pass
                        mykey = '{}_{}'.format(EntA,EntB)
                        if ComW > 0:
                            if mykey not in Com_dict:
                                Com_dict[mykey] = ComW
                            else:
                                prev_value = Com_dict[mykey]
                                if ComW > prev_value:
                                    Com_dict[mykey] = ComW
                                else:
                                    pass
                        if ExpW > 0:
                            if mykey not in Exp_dict:
                                Exp_dict[mykey] = ExpW
                            else:
                                prev_value = Exp_dict[mykey]
                                if ExpW > prev_value:
                                    Exp_dict[mykey] = ExpW
                                else:
                                    pass
            except KeyError:
                continue

    return Com_dict, Exp_dict

def make_edgelist(Com_dict,Exp_dict,File_Com,File_Exp): 
    for item in [(Com_dict,File_Com),(Exp_dict,File_Exp)]:
        mylist = [(key.split('_')[0],key.split('_')[1],float(val)/1000) for key, val in item[0].items()]
        mylist = sorted(mylist, key=lambda element: (element[0], element[1]))
        mydf = pd.DataFrame(mylist,columns=["Gene1","Gene2","Weight"])
        mydf.to_csv(item[1],sep="\t",header=True,index=False)
        
# set some filenames
fp_net = "../../data/networks/"
File_org = fp_net + '9606.protein.links.detailed.v12.0.txt'
File_Com = fp_net + 'STRING.edg'
File_Exp = fp_net + 'STRING-EXP.edg'
# run the functions
convert_dict = load_convert_dict()
Com_dict, Exp_dict = main(File_org)
make_edgelist(Com_dict,Exp_dict,File_Com,File_Exp)