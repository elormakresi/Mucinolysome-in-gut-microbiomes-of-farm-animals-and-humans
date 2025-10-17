import sys

##can ignore this script. I found a more efficient way to find genomes with coh/doc. This is only useful for the UHGG search
#python 3 selectbestdomain.py  cohesin.uhgp-100.domain.tbl dockerin_1.uhgp-100.domain.tbl genomes-all_metadata.tsv
## usage: python3 jerry_bestdomain_metegenome_information.py cohesin.uhgp-100.domain.tbl dockerin_1.uhgp-100.domain.tbl genomes-all_metadata.tsv

## python3 jerry_bestdomain_from_jerry.py cohesin.uhgp-100.domain.tbl dockerin_1.uhgp-100.domain.tbl genomes-all_metadata.tsv
## python3 jerry_bestdomain_from_jerry.py cohesin-new.uhgp-100.domain.tbl dockerin-new.uhgp-100.domain.tbl genomes-all_metadata.tsv


query2pos = {} ; query2cov = {}; query2domain = {}

# target name          accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
## GUT_GENOME269304_01280 -            482 Dockerin_1           PF00404.20    58   3.6e-60  212.8  53.2   1   5   6.8e-15   6.2e-12   58.3   1.0     2    58   109   164   108   164 0.96 hypothetical protein
import numpy as np
import pandas as pd
df = pd.read_csv("genomes-all_metadata.tsv",sep="\t",low_memory=False)
#df = pd.read_csv(sys.argv[3],sep="\t",low_memory=False)

for line in open(sys.argv[1]):
    if line.startswith("#"):
        continue
    lines = line.split()
    target_name = lines[0]
    start = lines[19]
    end  = lines[20]
    #print (target_name,start,end)
    cov = (int(end)-int(start) +1)/int(lines[2])
    query2pos.setdefault(target_name,[]).append(target_name+":"+start+":"+end)
    query2cov.setdefault(target_name,[]).append(cov)
    query2domain.setdefault(target_name,[]).append(lines[3])


query2pos_cohesin = {} ; query2cov_cohesin = {}; query2domain_cohesin = {}
for line in open(sys.argv[2]):
    if line.startswith("#"):
        continue
    lines = line.split()
    target_name = lines[0]
    start = lines[19]
    end  = lines[20]
    #print (target_name,start,end)
    cov = (int(end)-int(start) +1)/int(lines[2])
    query2pos_cohesin.setdefault(target_name,[]).append(target_name+":"+start+":"+end)
    query2cov_cohesin.setdefault(target_name,[]).append(cov)
    query2domain_cohesin.setdefault(target_name,[]).append(lines[3])

genome2proteinid = {} ; genome2domain = {}

#print ("Proteinid domain pos cov genomeid")

kept_genomeid = []

for target in query2cov:
    bestdomain_index =  np.argsort(query2cov[target])[-1]
    Proteinid,domain,pos,cov,genomeid = [target,query2domain[target][bestdomain_index],query2pos[target][bestdomain_index],query2cov[target][bestdomain_index],"_".join(target.split("_")[0:1])]
    #Proteinid,domain,pos,cov,genomeid = [target,query2domain[target][bestdomain_index],query2pos[target][bestdomain_index],query2cov[target][bestdomain_index],"_".join(target.split("_")[0:2])]
    
    
    genome2proteinid.setdefault(genomeid,[]).append(Proteinid)
    genome2domain.setdefault(genomeid,[]).append(domain)
for target in query2cov_cohesin:
    bestdomain_index =  np.argsort(query2cov_cohesin[target])[-1]
    Proteinid,domain,pos,cov,genomeid = [target,query2domain_cohesin[target][bestdomain_index],query2pos_cohesin[target][bestdomain_index],query2cov_cohesin[target][bestdomain_index],"_".join(target.split("_")[0:1])]
    genome2proteinid.setdefault(genomeid,[]).append(Proteinid)
    genome2domain.setdefault(genomeid,[]).append(domain)

for genomeid in genome2domain:
    #if "Cohesin" in genome2domain[genomeid] and "Dockerin_1" in genome2domain[genomeid]:
    Cohesin_number  = genome2domain[genomeid].count("Cohesin")
    Dockerin_number  = genome2domain[genomeid].count("Dockerin_1")
    
    if Cohesin_number >= 1 and Dockerin_number >= 1:
        #print (genomeid,genome2proteinid[genomeid])
        kept_genomeid.append(genomeid)
        values = df[df["Genome"] == genomeid].values
        if len(values) >= 1:
            value = values[0]
            print ("\t".join([str(ii) for ii in value]))
