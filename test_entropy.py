#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 16:01:04 2021

@author: yasmmin
"""
import pandas as pd
from datetime import datetime
import json

def run():
    sr={}
    with open('data/snap_data.json', 'r') as fp:
        sr = json.load(fp) 
        
    lineage="B.1.617.2"
    lineage="P.1"
    lineage="B.1.1.7"
    location="All"
    if(location=="All"):
        df=pd.read_csv('data/report_position_probabilities_by_lineage_period.tsv', sep='\t')
    #else:
    #    df=pd.read_csv(directory+"foca_backend/data/entropy/"+location.replace(" ", "-")+"_report_position_probabilities_by_lineage_period.tsv", sep='\t')
    df1 = df[ df["Lineage"] == lineage ]
    
    data={}
    labels=[]
    colors={}
    """
    for m in range(1,13):
        for y in [2020, 2021]:
            df=df[ ( (df['Month']==m) & (df['Year']==y) ) ]
            if(len(df)>0):
                df=df.sort_values('Entropy', ascending=False)
                
                y=df.iloc[0,5]
                if(y > 0):
                    xlab=str(m)+"/"+str(y)+" - "+str(df.iloc[0,4])
                    
                    prot=df.iloc[0,3]
                    if(not prot in data.keys()):
                        labels.append(prot)
                        
                        while (len(colors.values()) != len(labels)):
                            color=str(fake.hex_color())
                            if(not color in colors.values()):
                                colors[prot]=color
                                
                        data[prot]=[ [], [], { 'color': colors[prot] } ]
                    
                    data[prot][0].append(xlab)
                    data[prot][1].append(y)
                    """
     
    the=1.5
    mo=datetime.now().month
    y=datetime.now().year
    back=['#39C9CF','#3180D4', '#20b2aa', '#bc8f8f', "#E6E6FA", "#C68FDA"]
    co=0
    # check positions that maintain hight entropy
    posi={}
    ids=[]
    for m in range(mo-5, mo+1):
        df=df1[ ( (df1['Month']==m) & (df1['Year']==y) ) ]
        if(len(df)>0):
            lab=str(m)+"/"+str(y)
            posi[lab]={}
            df=df.sort_values('Entropy', ascending=False)  
            pos=[]              
            for i in range(len(df)):
                y_=df.iloc[i,8] # for count
                y_=df.iloc[i,5] # for entropy
                if(y_ > the):
                    p=str(df.iloc[i,4])
                    if(not p in pos):
                        pos.append(p)
                        prot=df.iloc[i,3]
                        posi[lab][prot+"_"+p]=y_
                        if(not prot+"_"+p in ids):
                            ids.append(prot+"_"+p)
    cnt={}
    for i in ids:
        cnt[i]=0
        for p in posi.keys():
            if(i in posi[p].keys()):
                cnt[i]+=1
        #if(cnt[i]>2):
        #    print(i,cnt[i])
    #print(cnt)
    
    eff={}
    neu={}
    co=0
    for m in range(mo-5, mo+1):
        df=df1[ ( (df1['Month']==m) & (df1['Year']==y) ) ]
        if(len(df)>0):
            df=df.sort_values('Entropy', ascending=False)
            lab=str(m)+"/"+str(y)
            
            """l
            if(not lab in data.keys()):
                labels.append(lab)
                
                while (len(colors.values()) != len(labels)):
                    color=str(fake.hex_color())
                    if(not color in colors.values()):
                        colors[lab]=color"""
                        
            data[lab]=[ [], [], { 'color': back[co] } ]
            
            eff[lab]={}
            neu[lab]={}
            
            pos = []
            end = len(df)
            if(end>20):
                end=20
                
            for i in range(len(df)):
                y_=df.iloc[i,8] # for count
                y_=df.iloc[i,5] # for entropy
                if(y_ > the):
                    p=str(df.iloc[i,4])
                    prot=df.iloc[i,3]
                    aas = df.iloc[i,7].split(" to ")
                    
                    sn=["-","-"]
                    if( prot.lower()+"_"+aas[0]+p+aas[1] in sr.keys()):
                        sn=sr[prot.lower()+"_"+aas[0]+p+aas[1]]
                        
                    xlab=prot+"_"+p
                    if(not p in pos):
                        if(xlab in cnt.keys()):
                            if(cnt[xlab]>2):
                                pos.append(p)
                                
                                data[lab][0].append(xlab)
                                data[lab][1].append(str(y_))
                                
                                eff[lab][xlab]=0
                                neu[lab][xlab]=0
                    
                    if(xlab in eff[lab].keys()):
                        if(sn[0] == "effect"):          
                            eff[lab][xlab]+=1
                        if(sn[0] == "neutral"):          
                            neu[lab][xlab]+=1
                
                if(len(pos)==6):
                    break
    print(data)
    print(eff)
    print(neu)
            
run()