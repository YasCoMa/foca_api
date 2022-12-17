#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 22:25:12 2021

@author: yasmmin
"""
directory="/var/www/html/foca_backend/"

import os
import json

import smtplib
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email import encoders

import urllib.request 

def _get_mutations(id_, seq, name_protein):
        f=open("temp.fasta","w")
        f.write(">"+str(id_)+"\n"+seq+"\n")
        f.close()
        
        mutations=""
        if(os.path.isfile(directory+"structural_effects/pipeline_strucomparison/blastdb/"+name_protein.lower()+"_ref.phr")):
            os.system("blastp -db "+directory+"structural_effects/pipeline_strucomparison/blastdb/"+name_protein.lower()+"_ref -query temp.fasta -num_alignments 1 -out temp.out -outfmt 3")
            
            refseq=""
            stref=0
            
            qseq=""
            stqseq=1
            
            found=False
            init=0
            co=0
            
            changes=[]
            f=open(directory+"temp.out","r")
            for line in f:
                if(line.startswith("Length")):
                    aux=line.replace("\n","").split("=")[1]
                    init=len(aux)+12
                    
                if(line.startswith("0 ") and not found):
                    l=line.replace("\n","").split(" ")
                    elements=[]
                    for el in l:
                        if(el!=""):
                            elements.append(el)
                    stref=int(elements[1])
                    refseq=line.replace("\n","")[init:init+60]
                    refseq=str(elements[2])
            
                    stqseq=0
                    for s in qseq:
                        if(s==' '):
                            break
                        else:
                            if(s!='X' and s!='N' and s!='n'):
                                if(refseq[stqseq] != "." and refseq[stqseq] != " " and (s!='-' or refseq[stqseq] != "-") ):
                                    changes.append( (refseq[stqseq])+str(stref)+s )
                            stref+=1
                            stqseq+=1
                            
                    found=True
                    
                if(line.startswith("Query_1")):
                    found=False
                    l=line.replace("\n","").split(" ")
                    elements=[]
                    for el in l:
                        if(el!=""):
                            elements.append(el)
                    qseq=str(elements[2])
                    #qseq=line.replace("\n","")[init:init+60]
                    
                co+=1
            f.close()
            
            mutations=';'.join(changes)
        return mutations
 
def _get_domain(sequence, protein):
    with open(directory+'data/train_domains.json', 'r') as fp:
        domains = json.load(fp)  
        
    doms=[]
    coords=[]
    if(protein in domains.keys()):
        for d in domains[protein]:
            c=0
            for s in domains[protein][d][1]:
                start=domains[protein][d][0][c][0]
                end=domains[protein][d][0][c][1]
                
                if(start<len(sequence)):
                    if(sequence[start-1:end]==s):
                        if(not d in doms):
                            doms.append(d)
                            coords.append(str(start)+"-"+str(end))
                c+=1
    
    c=0
    text=[]
    for d in doms:
        text.append(d+" ("+coords[c]+")")
        c+=1
        
    return ["; ".join(text), doms, coords]

def _get_muts_in_domains(sequence, coord, pos, muts):
    start = int(coord.split("-")[0])
    end = int(coord.split("-")[1])
    
    md=[]
    if(len(sequence) <= start):
        md.append('truncated_region')
        
    c=0
    for p in pos:
        if(int(p)>=start and int(p)<=end):
            md.append(muts[c])
    return md
    
def _search_domain(dat, id_):
    for d in dat.keys():
        if(dat[d]['prosite']==id_):
            return d
        
def get_domains_api(seq):
    url="https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi?seq=>seq%0A"+seq+"&output=txt"
    response = urllib.request.urlopen( url)
    doms = response.read() 
    
    with open(directory+'data/info_domains.json', 'r') as fp:
        dat = json.load(fp)
            
    text=[]
    domains=[]
    coords=[]
    doms=doms.splitlines()
    for line in doms:
        if(line!=""):
            line=str(line).replace("b'","").replace("'","")
            l=line.replace("\n","").split("\\t")
            if(len(l)>1):
                name=_search_domain(dat, l[3])
                
                text.append(name+" ("+l[1]+"-"+l[2]+")")
                domains.append(l[3]+"-"+name)
                coords.append(l[1]+"-"+l[2])
                
    return ["; ".join(text), domains, coords]

def send_success_email(dest, link, file) :
    name=file.split("/")[-1]
    
    fromaddr = "foca.app.lncc@gmail.com"
    toaddr = dest
    msg = MIMEMultipart()
    msg['From'] = fromaddr
    msg['To'] = toaddr
    
    msg['Subject'] = "Your Job in FOCA App finished"
    body = "<p>Hello, FOCA user. <br /> The files related to your job are attached in this e-mail</a> </p>"
    msg.attach(MIMEText(body, 'html'))
    
    part = MIMEBase("application", "octet-stream")
    part.set_payload(open(file + ".zip", "rb").read())
    encoders.encode_base64(part)
    part.add_header("Content-Disposition", "attachment; filename=\"%s.zip\"" % (name))
    msg.attach(part)
    
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login(fromaddr, "lncc2019")
    text = msg.as_string()
    server.sendmail(fromaddr, toaddr, text)
    server.quit()

    return None
