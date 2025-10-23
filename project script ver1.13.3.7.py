"""
Created on Thu Jul 30 16:22:28 2020

@author: s135271
"""


#!/usr/bin/env python3
import sys ## python -m pip install sys
import requests ## python -m pip install requests
import os ## python -m pip install os
import shutil ## python -m pip install shutil
import glob ## python -m pip install glob
from datetime import datetime ## python -m pip install datetime

# asks for substrate for transporter
print("What substrate would you like to find transporters for ?")
substrate=input().lower()
substrate=[substrate]

# as a minimum, the user needs to enter a substrate
if substrate==[]:
    sys.exit("WARNING!: As a minimum, you need to enter a substrate that you want to look for." + '\n' + "Write a substrate and then either enter an organisms name or press enter.")

# asks for organism to search for transporters in. 
print("Do you wish to search in a specific organism ?" + '\n' + "write the name of the organism. If you just press enter, it does for all.")
organism=input()
if organism=='':
    organism='all'

# set time for process initiation
times=datetime.now().strftime('%d-%m-%Y_%H;%M') # formate: day-month-year_hour;minutes

# function to translate organism name into NCBI taxi-ID
def orgtran(name):
    taxi_ID=None
    with open('names.dmp', 'r') as f: # open taxonomy dump file from NCBI
        for line in f:
            words = line.split("\t")
            if name.lower()==words[2].lower():
                taxi_ID=int(words[0])
                break ## we don't need to look through the rest of the file, if we found our hit
    return taxi_ID;

string_api_url = "https://string-db.org/api"
output_format = "tsv-no-header"
method = "get_string_ids"

# establishes the search terms for transporters
TRANSPORTERS = ["transporter",
                "pump",
                "receptor",
                "porin",
                "exporter",
                "antiporter",
                "porter",
                "carrier",
                "uniporter",
                "importer",
                "channel",
                "transport",
                "symporter",
                "pore",
                "nanotube",
                "TNT",
                "cell-permeable peptide",
                "CPP",
                "translocator"]

##
## Set parameters
##
params = {

    "identifiers" : "\r".join([str(substrate)]), # your substrate
    "species" : orgtran(organism), # species NCBI identifier 
    "limit" : 100, # only the 100 (best) identifier per substrate
    "echo_query" : 1, # see your input identifiers in the output
    #"caller_identity" : "www.awesome_app.org" # your app name

}

##
## Construct URL
##
request_url = "/".join([string_api_url, output_format, method])


##
## Protein sequence-calling function
##
def protseq(string_ID):
    ## Make search in order to get accession ID, annotated name, 
    ##Pfam numbers and protein sequence in fasta formate
    url = 'https://www.uniprot.org/uploadlists/'
    params = {'from': 'STRING_ID', # search based on STRING ID
        'to': 'ACC', # search for UniProtKB accessions and find entry
        'format': 'txt', # get entry in text format
        'include': 'yes', # includes isoforms of the protein
        'query': string_ID # query for the STRING ID of the given entry
        }
    ## make output workable
    tex = requests.post(url,data=params)
    tjek0=tex.text.split('\n')
    
    local="NA"
    Pfam=[]
    FuncKey=[]
    seq=''
    flag="down"
    HITS=0
    flag1="down"
    TMBcheck="none"
    TMB=0
    FullName=""
    nameflag="access granted"
    annot=["TRANSMEM","TOPO_DOM","note","evidence"]
    if len(tjek0)>1: ## if there ARE results
        for line in tjek0:
            LINE=line.split()
            if len(LINE)>1: ## the last line is empty with len=1
                if LINE=='//': ## end sequence fetching
                    flag="down"
                if LINE[0]=='AC': ## extracting accession ID
                    acc=LINE[-1].strip(';')
                    header='>' + LINE[-1] + '\n' ## using accession ID as header
                    HITS+=1
                elif LINE[0]=='DE' or LINE[0]=='KW': ## UniProt annotation for match checking
                    FuncKey.extend(LINE)
                    if LINE[0]=='DE' and len(LINE)>2 and nameflag=="access granted":
                        if 'Full=' in LINE[2]:
                            del LINE[:2]
                            del LINE[-1]
                            FullNamet=' '.join(LINE)
                            FullNamet=FullNamet[5:]
                            for x in FullNamet:
                                if x!="{" and x!=";" and x!=",":
                                    FullName+=x
                                else:
                                    break
                            nameflag="access denied"
                elif LINE[0]=="FT":
                    if flag1=="up":
                        if TMBcheck=="check":
                            if not any(x in LINE[1] for x in annot):
                                TMB+=1
                                TMBcheck="none"
                        elif "Beta" in LINE[1]:
                            TMBcheck="check"
                    elif LINE[1]=="TRANSMEM":
                        flag1="up"
                elif LINE[0]=='CC': 
                    if 'LOCATION:' in LINE: # checks for subcellular location
                        local=""
                        del LINE[:4] ## only focuses on the informational part of the line
                        #print(LINE)
                        for x in ' '.join(LINE):
                            if x!="{" and x!=";": ## only saves the necessary information
                                local+=x
                            else:
                                break ## omae wa mou shinderu
                elif LINE[1]=='Pfam;': ## extracting Pfam domain codes
                    Pfam.append(LINE[2].strip(';'))
                if flag=="up": ## extracting sequence as fasta
                    linje=''.join(LINE) + '\n'
                    seq+=linje
                elif LINE[0]=='SQ': ## begin sequence extraction with next line
                    flag="up"
    else: ## if there are NO results
        acc="NA"
        FuncKey=[]
        header=None
        seq=None
        Pfam=["NA"]
        local="NA"
        TMB="NaN"
        FullName="NA"
    if header!=None and seq!=None: ## makes a fitting fasta-format
        entry=header+seq
        ## this format is easier to save to file for usage in prediction.
        FuncKey= [item.lower() for item in FuncKey] ## saves annotation as list to analyse later.
    else:
        entry=''
        
    return (entry,HITS,acc,FuncKey,Pfam,local,TMB,FullName);


##
## TMHMM function for TM-helox predictions
##
def tmhmm(fasta):
    ## write found fasta into .fasta
    ## since TMHMM requires it in order to run
    with open('test.fasta','w') as test:
        for line in fasta:
            test.write(line)
        test.close()

    ## run TMHMM and get output as output-variable "record"
    tmhmm='tmhmm-2.0c/bin/tmhmm ' + 'test.fasta'
    record=os.popen(tmhmm).read().split('\n')
    
    ## no more need for the written fasta-file, so it is deleted
    os.remove('test.fasta')
    
    ## TMHMM also generates a folder that is called TMHMM_*****
    ## we just wants the printed output, no the documents, so the folder is deleted
    TMTRASH=glob.glob('TMHMM_*')
    shutil.rmtree(''.join(TMTRASH))
    
    ## get the number og predicted TM-helices
    # start with empty list in order to get predictions for all output proteins
    TMHMMpred=[]
    for line in record:
        
        ## only do stuff, if the line actually tells us the number of predicted TM-helices
        if "Number of predicted TMHs" in line:
            ## In that line, the number is the last element in the line, 
            ##so splitting it by whitespace and calling the final element 
            ##in the line will call the number og predicted TM-helices
            LINE=line.split()
            TMHMMpred.append(LINE[-1])
            
    ## if there were no protein hits, it will give the value 'Nan'
    if TMHMMpred==[]:
        TMHMMpred='NaN'
    ## else if there are only 1 entry found, there are no reason to print the number as a list.  
    elif len(TMHMMpred)==1:
        TMHMMpred=''.join(TMHMMpred)
    
    return(TMHMMpred);

##
## Call STRING
##
results = requests.post(request_url, data=params)
nrHITS=0
nrUNI=0
spec=[]
TMHcount=0
matchcount=0
mismatchcount=0
animation="\|/\|/"
time=0
##
## Read and parse the results
##

outputTable=[["INPUT","ORGANISM","PROTEIN","#UP","#TMH","#TMBB","Location","type","accession ID","Protein Name","Pfam domains"]]

if  results.text!='': # only go through the analytical steps, if STRING actually can FIND ANYTHING
    for line in results.text.strip().split("\n"):
        l = line.split("\t")
        input_identifier, STRING, string_identifier, prot_identifier = ''.join(substrate), l[2], l[4], l[5]
        ## loading animation
        time+=1
        sys.stdout.write("\r" + "="+ "=" + animation[time % len(animation)] +
                         animation[time % len(animation)] + "=" + "=" +
                     "" + ' progress ' + str(int(time/len(results.text.strip().split("\n"))*100)) + "%" + " " +
                     "=" + "=" + animation[time % len(animation)] +
                     animation[time % len(animation)] + "=" + "=")
        sys.stdout.flush()
        
        if any(x in l[6].lower() for x in TRANSPORTERS): ## checks the STRING-annotation, if the entry is a transporter
            for x in TRANSPORTERS: ## goes through the annotation again in order to save the term checked 
                if x in l[6].lower():
                    transkey=x
                    break
            result,prothit,accnr,UNI,PfamID,place,TMBB,protName = protseq(STRING) ## acquires information from UniProt
            ## check to determine, if there are a match between the annotation of STRING and UniProt
            termcheck=""
            for x in UNI:
                if transkey in UNI:
                    termcheck="match"
            if termcheck=="":
                if prothit==0:    
                    termcheck="NA"
                else:
                    termcheck="mismatch"
            ## performs TMHMM
            TMH=tmhmm(result)
            outputTable.append([input_identifier, string_identifier, prot_identifier, str(prothit), str(TMH), str(TMBB), place, transkey, accnr, protName, ','.join(PfamID)])
            ## gathers data on the output like amount of STRING- and UniProt-hits and how many had TMH's
            nrHITS+=1
            if prothit>0:
                nrUNI+=1
            if string_identifier not in spec:
                spec.append(string_identifier)
            if prothit==0:
                test=l
            if termcheck=="match":
                matchcount+=1
            if termcheck=="mismatch":
                mismatchcount+=1
            if TMH!='NaN':
                if int(TMH)>0:
                    TMHcount+=1
    a_string = "{:<15} {:<40} {:<13} {:<3} {:<4} {:<5} {:<20} {:<12} {:<12} {:<54} {:<17}"
    
    print('\n')
    n=0
    # the output files will (for orders sake) be saved in an output folder. 
    # but since this folder is specific for these output files, it will check if 
    # the output folder exists already and make one, if it does not. 
    if not os.path.exists("output"):
        os.makedir("output")
    # set name for output file
    name="output/"+input_identifier+"_"+organism+"_"+str(times)+".txt" # format: substrate_organism_day-month-year_hour;minute
    with open(name,"w") as file:
        for line in outputTable:
            file.write('\t'.join(line)+'\n')
            print(a_string.format(*outputTable[n]))
            n+=1
    
    print('\n')
    # # prints the statistics on the output 
    print("{:<15} {:<15} {:<10} {:<25} {:<25} {:<30}".format("protein hits",
                                                             "UniProt hits",
                                                             "# species",
                                                             "# annotation matches",
                                                             "# annotation mismatches",
                                                             "# entries with pred TMH"))
    print("{:<15} {:<15} {:<10} {:<25} {:<25} {:<30}".format(nrHITS,
                                                         nrUNI,
                                                         len(spec),
                                                         matchcount,
                                                         mismatchcount,
                                                         TMHcount))
else: ## if no STIRNG-entries were found. 
    print("No STRING hits were found")

