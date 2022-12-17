#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 13:08:43 2021

@author: yasmmin
"""

import pandas as pd


from datetime import datetime
from datetime import date

import operator
import statistics
import math
import json
import os

from matplotlib import pyplot as plt
import fbprophet

class Generate_processed_data:
    
    def get_all_snap_results(self):
        snap={}
        folder="data/snap_results"
        for f in os.listdir("data/snap_results"):
            protein=f.replace(".csv","")
            
            c=0
            g=open(folder+"/"+f,"r")
            for line in g:
                if(c>0):
                    l=line.replace("\n","").split(",")
                    snap[protein+"_"+l[0]] = [l[1], l[2]]
                c+=1
            g.close()
        
        with open('data/snap_data.json', 'w') as fp:
            json.dump(snap, fp)
        
    def get_demography_data(self):
        locs=[]
        f=open("data/location.tsv","r")
        for line in f:
            locs.append(line.replace("\n",""))
        f.close()
        
        with open('data/status_epi.json', 'r') as fp:
            dat = json.load(fp)  
        
        names_vocs=["Gamma","Alpha","Beta","Delta","Omicron"]
        vocs=["P.1","B.1.1.7","B.1.351","B.1.617.2","B.1.1.529"]
        condsvc=["( lineage=='P.1' or lineage.startswith('P.1.') )", "lineage=='B.1.1.7'", "lineage=='B.1.351'", "( lineage=='B.1.617.2' or lineage.startswith('AY.') )"]
        condsvc=["( lineage=='P.1' or lineage.startswith('P.1.') or variant.find('VOC Gamma')!=-1 )", "(lineage=='B.1.1.7' or variant.find('VOC Alfa')!=-1 )", "(lineage=='B.1.351' or variant.find('VOC Beta')!=-1)", "( lineage=='B.1.617.2' or lineage.startswith('AY.') or variant.find('VOC Delta')!=-1 )", "( lineage=='B.1.1.529' or variant.find('VOC Omicron')!=-1 )"]
        
        conds=['( status.find("mild")!=-1 or status.find("live")!=-1 or status.find("home")!=-1 or status.find("not hosp")!=-1 or status.find("cured")!=-1 or status.find("outpatient")!=-1 or status.find("released")!=-1 or status.find("home")!=-1 )', '( status.find("asymptomatic")!=-1 )','( status.find("moderate")!=-1 )', '( status.find("hosp")!=-1 and status.find("not hosp")==-1 )', '( status.find("severe")!=-1 or status.find("severa")!=-1 or status.find("critical")!=-1 )', '( status.find("icu")!=-1 or status.find("uci")!=-1 or status.find("intensive")!=-1 )', '( status.find("death")!=-1 or status.find("deceased")!=-1 or status.find("died") or status.find("fatal")!=-1 or status.find("dead")!=-1 )']
        names_status=["mild/alive", "asymptomatic", "moderate", "hospitalized", "severe", "icu", "deceased"]
        
        conds=['( status.find("mild")!=-1 or status.find("live")!=-1 or status.find("not hosp")!=-1 or status.find("cured")!=-1 or status.find("outpatient")!=-1 or status.find("released")!=-1 or status.find("home")!=-1 )', '( status.find("hosp")!=-1 and status.find("not hosp")==-1 ) or ( status.find("severe")!=-1 or status.find("severa")!=-1 or status.find("critical")!=-1 ) or ( status.find("icu")!=-1 or status.find("uci")!=-1 or status.find("intensive")!=-1 )', '( status.find("death")!=-1 or status.find("deceased")!=-1 or status.find("died")!=-1 or status.find("dead")!=-1 )']
        names_status=["Mild/Alive", "Hospitalized", "Deceased"]
        
        tconds=['( status.find("mild")!=-1 or status.find("live")!=-1 or status.find("release")!=-1 or status.find("not hosp")!=-1 or status.find("cured")!=-1 or status.find("outpatient")!=-1 or status.find("released")!=-1 or status.find("home")!=-1 )', '( status.find("death")!=-1 or status.find("deceased")!=-1 or status.find("died")!=-1 or status.find("dead")!=-1 or status.find("severe")!=-1 or status.find("severa")!=-1 or status.find("critical")!=-1 or status.find("icu")!=-1 or status.find("uci")!=-1 or status.find("intensive")!=-1 )']
        tnames_status=["mild", "severe"]
        
        list_status=set()
        
        genders=["Male","Female"]
        ages=['Child/Adolescent','Adult','Elderly']
        filters=['samples', 'Child/Adolescent','Adult','Elderly']+["Male","Female"]+names_status+names_vocs
        global_summary={}
        for l in locs:
            global_summary[l]={}
            for fi in filters:
                global_summary[l][fi]=0
        
        age_by_voc = {}
        for v in ['Child/Adolescent','Adult','Elderly']:
            age_by_voc[v]={}
            for u in names_vocs:
                age_by_voc[v][u]=0
                
        age_by_status = {}
        for v in ['Child/Adolescent','Adult','Elderly']:
            age_by_status[v]={}
            for u in names_status:
                age_by_status[v][u]=0
                
        status_by_voc = {}
        for v in names_status:
            status_by_voc[v]={}
            for u in names_vocs:
                status_by_voc[v][u]=0
                
        by_age={"Male": {}, "Female": {}}
        for v in ['Child/Adolescent','Adult','Elderly']:
            by_age["Male"][v]=0
            by_age["Female"][v]=0
            
        by_voc={"Male": {}, "Female": {}}
        for v in names_vocs:
            by_voc["Male"][v]=0
            by_voc["Female"][v]=0
            
        by_patient_status={"Male": {}, "Female": {}}
        for v in names_status:
            by_patient_status["Male"][v]=0
            by_patient_status["Female"][v]=0
        
        states = {
    'AC': 'Acre',
    'AL': 'Alagoas',
    'AP': 'Amapa',
    'AP_': 'Amapá',
    'AM': 'Amazonas',
    'BA': 'Bahia',
    'CE': 'Ceara',
    'CE_': 'Ceará',
    'DF': 'Distrito Federal',
    'ES': 'Espirito Santo',
    'ES_': 'Espírito Santo',
    'GO': 'Goias',
    'GO_': 'Goiás',
    'MA': 'Maranhao',
    'MA_': 'Maranhão',
    'MT': 'Mato Grosso',
    'MS': 'Mato Grosso do Sul',
    'MG': 'Minas Gerais',
    'PA': 'Para',
    'PA_': 'Pará',
    'PB': 'Paraiba',
    'PB_': 'Paraíba',
    'PR': 'Parana',
    'PR_': 'Paraná',
    'PE': 'Pernambuco',
    'PI': 'Piaui',
    'PI_': 'Piauí',
    'RJ': 'Rio de Janeiro',
    'RN': 'Rio Grande do Norte',
    'RS': 'Rio Grande do Sul',
    'RO': 'Rondonia',
    'RO_': 'Rondônia',
    'RR': 'Roraima',
    'SC': 'Santa Catarina',
    'SP': 'Sao Paulo',
    'SP_': 'São Paulo',
    'SE': 'Sergipe',
    'TO': 'Tocantins'
}
        mapst={}
        for k in states:
            mapst[states[k].lower()]=k.replace("_","")
        mapstrev={}
        
        for k in states:
            mapstrev[k.replace("_","")]=states[k]
            
        state_mutations = {}
        state_lineages = {}
        count_state = {}
        
        g=open("data/demography/metadata_demography.tsv","w")
        g.write("VOC\tEPI\tLocation\tLineage\tGender\tAge\tStatus\tCollection_Date\n")
        
        c=0
        f=open("data/metadata.tsv","r")
        for line in f:
            if(c>0):
                l=line.replace("\n","").split("\t")
                epi=l[2]
                dt=l[3]
                lineage=l[11]
                variant = l[13]
                gender=l[9]
                voc=""
                
                try:
                    age=int(l[8])
                except:
                    try:
                        a1=int(l[8].split("-")[0])
                        a2=int(l[8].split("-")[1])
                        age=int( (a1+a2)/2 )
                    except:
                        age=-1
                        
                location=""
                for lc in locs:
                    if(l[4].find(lc)!=-1):
                        location=lc
                        break
                
                status="unknown"
                if(epi in dat.keys()):
                    status=dat[epi].lower()
                    
                    if(not dat[epi] in list_status):
                        list_status.add(dat[epi].lower())
                
                if(location!=""):
                    count=0
                    for cd in condsvc:
                        if(eval(cd)):
                            voc=names_vocs[count]
                        count+=1
                        
                    if(status!="unknown"):
                        global_summary[location]['samples']+=1
                    
                        if(age>0 and age<18):
                            global_summary[location]['Child/Adolescent']+=1
                        if(age>=18 and age<=60):
                            global_summary[location]['Adult']+=1
                        if(age>60):
                            global_summary[location]['Elderly']+=1
                            
                        count=0
                        for cd in conds:
                            if(eval(cd)):
                                global_summary[location][names_status[count]]+=1
                            count+=1
                            
                        count=0
                        for cd in condsvc:
                            if(eval(cd)):
                                global_summary[location][names_vocs[count]]+=1
                            count+=1
                            
                        #if(lineage in vocs):
                        #    global_summary[location][lineage]+=1
                        
                        if(gender in ['Female','Male']):
                            global_summary[location][gender]+=1
                        
                if(location=="Brazil"):
                    if(len(l[4].split(" / "))>2):
                        if(l[4].split(" / ")[2].lower() in mapst.keys()):
                            state = "BR-"+mapst[l[4].split(" / ")[2].lower()]
                            if(not state in state_mutations.keys()):
                                state_lineages[state] = {}
                                state_mutations[state] = {}
                                count_state[state]=0
                            count_state[state]+=1
                            
                            if(lineage!="unknown" and lineage!=""):
                                if(not lineage in state_lineages[state].keys()):
                                    state_lineages[state][lineage]=0
                                state_lineages[state][lineage]+=1
                            try:
                                if(l[14]!="" and l[14]!="NA"):
                                    mutations=l[14].replace("(","").replace(")","").split(",")
                                    for m in mutations:
                                        if(not m in state_mutations[state].keys()):
                                            state_mutations[state][m]=0
                                        state_mutations[state][m]+=1
                            except:
                                pass
                    
                    # Creating status_by_voc
                    if(status != "unknown"):
                        count=0
                        for cd in conds:
                            if(eval(cd)):
                                 count2=0
                                 for cd2 in condsvc:
                                    if(eval(cd2)):
                                        status_by_voc[names_status[count]][names_vocs[count2]]+=1
                                    count2+=1
                                    
                                #if(lineage in vocs):
                                #    status_by_voc[names_status[count]][lineage]+=1
                            count+=1
                    
                    # Creating age_by_voc and age_by_status
                    if(age>0 and age<18):
                        count=0
                        for cd in condsvc:
                            if(eval(cd)):
                                age_by_voc['Child/Adolescent'][names_vocs[count]]+=1
                            count+=1
                            
                        #if(lineage in vocs):
                        #    age_by_voc['Child_Adolescent'][lineage]+=1
                        
                        if(status != "unknown"):
                            count=0
                            for cd in conds:
                                if(eval(cd)):
                                    age_by_status['Child/Adolescent'][names_status[count]]+=1
                                count+=1
                            
                    if(age>=18 and age<=60):
                        count=0
                        for cd in condsvc:
                            if(eval(cd)):
                                age_by_voc['Adult'][names_vocs[count]]+=1
                            count+=1
                            
                        #if(lineage in vocs):
                        #    age_by_voc['Adult'][lineage]+=1
                        
                        if(status != "unknown"):
                            count=0
                            for cd in conds:
                                if(eval(cd)):
                                    age_by_status['Adult'][names_status[count]]+=1
                                count+=1
                                
                    if(age>60):
                        count=0
                        for cd in condsvc:
                            if(eval(cd)):
                                age_by_voc['Elderly'][names_vocs[count]]+=1
                            count+=1
                            
                        #if(lineage in vocs):
                        #    age_by_voc['Elderly'][lineage]+=1
                        
                        if(status != "unknown"):
                            count=0
                            for cd in conds:
                                if(eval(cd)):
                                    age_by_status['Elderly'][names_status[count]]+=1
                                count+=1
                                
                    # Creating by_age, by_voc and by_patient_status    
                    if(gender in ['Female','Male']):
                        if(age>0 and age<18):
                            by_age[gender]['Child/Adolescent']+=1
                        if(age>=18 and age<=60):
                            by_age[gender]['Adult']+=1
                        if(age>60):
                            by_age[gender]['Elderly']+=1
                        
                        
                        count=0
                        for cd in condsvc:
                            if(eval(cd)):
                                by_voc[gender][names_vocs[count]]+=1
                            count+=1
                            
                        #if(lineage in vocs):
                        #    by_voc[gender][lineage]+=1
                            
                        if(status != "unknown"):
                            count=0
                            for cd in conds:
                                if(eval(cd)):
                                    by_patient_status[gender][names_status[count]]+=1
                                count+=1
                                
                            #if(not status in by_patient_status[gender].keys()):
                            #    by_patient_status[gender][status]=0
                            #by_patient_status[gender][status]+=1
                        
                if(status != "unknown"):
                    st=[]
                    count=0
                    for cd in tconds:
                        if(eval(cd)):
                            st.append(tnames_status[count])
                        count+=1
                    status = ';'.join(st)
                    
                g.write("%s\t%s\t%s\t%s\t%s\t%i\t%s\t%s\n" %(voc, epi, location, lineage, gender, age, status, dt) )
            c+=1
        f.close()
        g.close()
        
        global_data=[]
        global_data.append(["Country", "Number of Samples"])
        global_info={}
        mapp={}
        df = pd.read_csv('data/codes.csv')
        for i in range(len(df)):
            co=df.iloc[i,0].split(",")[0]
            if(df.iloc[i,0]=="Virgin Islands, British"):
                co="British Virgin Islands"
            if(df.iloc[i,0]=="Virgin Islands, U.S."):
                co="Virgin Islands"
            if(df.iloc[i,0]=="Brunei Darussalam"):
                co="Brunei"
            if(df.iloc[i,0]=="Korea, Republic of"):
                co="South Korea"
            if(df.iloc[i,0]=="Cape Verde"):
                co="Cabo Verde"
            if(df.iloc[i,0]=="Côte d'Ivoire"):
                co="Cote d'Ivoire"
            if(df.iloc[i,0]=="Curaçao"):
                co="Curacao"
            if(df.iloc[i,0]=="Congo, the Democratic Republic of the"):
                co="Democratic Republic of the Congo"
            if(df.iloc[i,0]=="Guinea-Bissau"):
                co="Guinea Bissau"
            if(df.iloc[i,0]=="Libya"):
                co="Libia"
            if(df.iloc[i,0]=="Russian Federation"):
                co="Russia"
            if(df.iloc[i,0]=="Saint Barthélemy"):
                co="Saint Barthelemy"
            if(df.iloc[i,0]=="Saint Martin (French part)"):
                co="Saint Martin"
            if(df.iloc[i,0]=="Sint Maarten (Dutch part)"):
                co="Sint Maarten"
            if(df.iloc[i,0]=="Viet Nam"):
                co="Vietnam"
            
            mapp[co]=str(df.iloc[i,1])
        
        with open("data/demography/region_countries.json", 'w') as fp:
            json.dump(mapp, fp)
    
        for k in global_summary.keys():
            cou=k
            if(k=="USA"):
                cou="United States"
            if(cou in mapp.keys()):
                if(global_summary[k]['samples']>0):
                    sub=[cou, global_summary[k]['samples']]
                    global_data.append( sub )
                    
                    #global_info[mapp[cou]]="<b>Infomations about "+cou+":</b><br /> <ul>" 
                    global_info[mapp[cou]]='''
                    <div style="border: 2px solid #000; text-align: center; border-radius: 10px;">
                    <p style="font-size: 20px; padding: 15px;"> Informations about '''+cou+'''</p>
                    '''
                    
                    ndimensions=['Age', 'Gender', 'Patient Status', 'VOC']
                    dimensions=[ages, genders, names_status, names_vocs]
                    
                    c=0
                    for d in dimensions:
                        global_info[mapp[cou]]+='<p style="border-top :1px solid #000; border-bottom :1px solid #000; padding: 5px 10px;"> By '+ndimensions[c]+': </p>'
                        temp_sum=0
                        for v in d:
                            vr=(global_summary[k][v] / global_summary[k]['samples'])*100
                            global_info[mapp[cou]]+='<div style="padding: 10px"><b>'+(v+":</b> "+str(global_summary[k][v])+f' ({vr:.2f}%) </div>')
                            temp_sum+=global_summary[k][v]
                        
                        vr=( (global_summary[k]['samples'] - temp_sum) / global_summary[k]['samples'])*100
                        global_info[mapp[cou]]+='<div style="padding: 10px"><b>'+("Other:</b> "+str(temp_sum)+f' ({vr:.2f}%) </div>')
                        c+=1
                        
                    """c=1
                    for v in list(global_summary[k].values())[1:]:
                        vr=(v/global_summary[k]['samples'])*100
                        global_info[mapp[cou]]+="<li> "+(filters[c]+": "+str(v)+f' ({vr:.2f}%) </li>')
                        c+=1"""
                        
                    global_info[mapp[cou]]+="</div>" 
            else:
                print(cou)
        with open("data/demography/global_data.json", 'w') as fp:
            json.dump(global_data, fp)
        with open("data/demography/global_info.json", 'w') as fp:
            json.dump(global_info, fp)
            
        f=open("data/demography/list_status_all.tsv","w")
        for l in list_status:
            f.write(l+"\n")
            
            """if(not l in by_patient_status['Male'].keys()):
                by_patient_status['Male'][l]=0
            if(not l in by_patient_status['Female'].keys()):
                by_patient_status['Female'][l]=0"""
        f.close()
        
        f=open("data/demography/list_status.tsv","w")
        for l in names_status :
            f.write(l+"\n")
        f.close()
        
        back=['#DD8673','#73ADDD', '#20b2aa', '#bc8f8f', "#E6E6FA"]
        back=['0A2239','53A2BE','1D84B5','901483','C76B99','982649','FB4B5A','D08C60','E05200','849263','FCD422']
        border=['lightred','lightblue', '#008b8b', '#a0522d','#D8BFD8']
        names=['by_age', 'by_voc', 'by_patient_status', 'age_by_voc', 'age_by_status', 'status_by_voc']
        dicts=[by_age, by_voc, by_patient_status, age_by_voc, age_by_status, status_by_voc]
        co=0
        for dc in dicts:
            first_key = list(dc.keys())[0]
            data_={"labels": list(dc[first_key].keys()), "datasets": [] }
            c=0
            for k in dc.keys():
                ds={}
                ds['label']=k
                ds['backgroundColor']="#"+back[c]
                ds['borderColor']="#"+back[c]
                ds['borderWidth']=1
                ds['data']=list(dc[k].values())
                
                data_['datasets'].append(ds)
                c+=1
                
            with open("data/demography/"+names[co]+'_data.json', 'w') as fp:
                json.dump(data_, fp)
                
            co+=1
        
        
        back=['#DD8673','#73ADDD', '#20b2aa', '#bc8f8f', "#E6E6FA"]
        back=['0A2239','53A2BE','1D84B5','901483','C76B99','982649','FB4B5A','D08C60','E05200','849263','FCD422']
        border=['#1D5557','#193B5F', '#008b8b', '#a0522d','#D8BFD8']
        names=[ 'state_lineages', 'state_mutations' ]
        dicts=[state_lineages , state_mutations ]
        co=0
        dss={}
        for dc in dicts:
            dss[names[co]]={}
            first_key = list(dc.keys())[0]
            
            for k in dc.keys():
                vs=sorted(dc[k].items(),key=operator.itemgetter(1),reverse=True)
                labels=[]
                values=[]
                for kv in vs:
                    if(len(labels)<5):
                        labels.append(kv[0])
                        values.append(kv[1])
                        
                data_={"labels": labels, "datasets": [] }
                
                ds={}
                ds['label']=k
                ds['backgroundColor']="#"+back[co]
                ds['borderColor']="#"+back[co]
                ds['borderWidth']=1
                ds['data']=values
                
                data_['datasets'].append(ds)
                
                dss[names[co]][k]=data_
            co+=1
        
        geostatedata=[["State","Samples"]]
        for k in count_state.keys():
            geostatedata.append( [{"v": k, "f": mapstrev[k.replace("BR-","")] }, count_state[k]] )
            
        with open("data/br_state_analysis_plot.json", 'w') as fp:
            json.dump(dss, fp)
        with open("data/br_state_analysis_data.json", 'w') as fp:
            json.dump(geostatedata, fp)
        
    def get_mutation_indel_data_lineage(self):
        locs=[]
        f=open("data/lineage.tsv","r")
        for line in f:
            locs.append(line.replace("\n",""))
        f.close()
        
        proteins=[]
        f=open("data/proteins.tsv","r")
        for line in f:
            proteins.append(line.replace("\n",""))
        f.close()
        
        df=pd.read_csv("data/metadata.tsv", sep="\t")
        names=["weekly","monthly"]
        co=0
        for period in ["w%V/%Y", "%m/%Y"]:
            data={"indel": {}, "missense": {}}
            
            for l in locs:
                #if(not os.path.isfile('data/lineage_mutation_indel/'+l.replace(" ","-")+"_"+names[co]+'_mutation_indel.json')):
                aux = df[ df['Pango lineage']==l ]
                aux['Collection date'] = pd.to_datetime(aux["Collection date"])
                aux=aux.sort_values(by='Collection date')
                ant=""
                labels=[]
                
                indel={}
                mis={}
                for p in proteins:
                    indel[p]=[]
                    mis[p]=[]
                    
                    data['missense'][p]={"data": [], "deviation": []}
                    data['indel'][p]={"data": [], "deviation": []}
                    
                for i in range(len(aux)):
                    mutsa={}
                    mutsb={}
                    for p in proteins:
                        mutsa[p]=[]
                        mutsb[p]=[]
                    
                    try:
                        if( len(str(aux.iloc[i,3]).split("-")) > 2 ):
                            #dat = datetime.fromisoformat(aux.iloc[i,3])
                            dat=aux.iloc[i,3]
                            if( int(dat.strftime("%Y")) > 2019):
                                
                                if(aux.iloc[i,14]!="NA" and aux.iloc[i,14]!=""):
                                    try:
                                        mutations=aux.iloc[i,14].replace("(","").replace(")","").split(",")
                                   
                                        for m in mutations:
                                        
                                            pr=m.split("_")[0]
                                            m=m.split("_")[1]
                                            if(m.find("del")!=-1 or m.find("ins")!=-1):
                                                mutsa[pr].append(m)
                                            else:
                                                mutsb[pr].append(m)
                                    except:
                                        pass
                                        
                                    for pr in proteins:
                                        indel[pr].append( len(mutsa[pr]) )
                                        mis[pr].append( len(mutsb[pr]) )
                                else:
                                    for pr in proteins:
                                        indel[pr].append(0)
                                        mis[pr].append(0)
            
                                per=dat.strftime(period)
            
                                if(ant!=per):
                                    if(ant!=""):
                                        labels.append(ant)
                                        for p in proteins:
                                            if( len(indel[p])==0 ):
                                                indel[p].append(0)
                                            if( len(indel[p])==1 ):
                                                indel[p].append(0)
                                                
                                            mean = statistics.median(indel[p])
                                            st = statistics.stdev(indel[p])
                                            data['indel'][p]["data"].append( mean )
                                            data['indel'][p]["deviation"].append( st )
                                            
                                            if( len(mis[p])==0 ):
                                                mis[p].append(0)
                                            if( len(mis[p])==1 ):
                                                mis[p].append(0)
                                                
                                            mean = statistics.median(mis[p])
                                            st = statistics.stdev(mis[p])
                                            data['missense'][p]["data"].append( mean )
                                            data['missense'][p]["deviation"].append( st )
                                    ant=per
                                    for p in proteins:
                                        indel[p]=[]
                                        mis[p]=[]
            
                    except:
                       pass
                    
                sum_=0
                for p in proteins:
                    sum_+=len(indel[p])
                
                if(sum_ > 0):
                    labels.append(ant)
                    for p in proteins:
                        if( len(indel[p])==0 ):
                            indel[p].append(0)
                        if( len(indel[p])==1 ):
                            indel[p].append(0)
                            
                        mean = statistics.median(indel[p])
                        st = statistics.stdev(indel[p])
                        data['indel'][p]["data"].append( mean )
                        data['indel'][p]["deviation"].append( st )
                        
                sum_=0
                for p in proteins:
                    sum_+=len(mis[p])
                
                if(sum_ > 0):
                    for p in proteins:
                        if( len(mis[p])==0 ):
                            mis[p].append(0)
                        if( len(mis[p])==1 ):
                            mis[p].append(0)
                            
                        mean = statistics.median(mis[p])
                        st = statistics.stdev(mis[p])
                        data['missense'][p]["data"].append( mean )
                        data['missense'][p]["deviation"].append( st )
            
                data['labels'] = labels 
                with open('data/lineage_mutation_indel/'+l.replace(".","-")+"_"+names[co]+'_mutation_indel.json', 'w') as fp:
                    json.dump(data, fp)
            
            co+=1
    
    def get_mutation_indel_data_general_location(self):
        locs=[]
        f=open("data/location.tsv","r")
        for line in f:
            locs.append(line.replace("\n",""))
        f.close()
        
        proteins=[]
        f=open("data/proteins.tsv","r")
        for line in f:
            proteins.append(line.replace("\n",""))
        f.close()
        
        df=pd.read_csv("data/metadata.tsv", sep="\t")
        names=["weekly","monthly"]
        co=0
        for period in ["w%V/%Y", "%m/%Y"]:
            data={"indel": {}, "missense": {}}
            
            for l in locs:
                #if(not os.path.isfile('data/mutation_indel/'+l.replace(" ","-")+"_"+names[co]+'_mutation_indel.json')):
                aux = df[ df['Location'].str.contains(l, na=False) ]
                aux['Collection date'] = pd.to_datetime(aux["Collection date"])
                aux=aux.sort_values(by='Collection date')
                ant=""
                labels=[]
                
                indel={}
                mis={}
                for p in proteins:
                    indel[p]=[]
                    mis[p]=[]
                    
                    data['missense'][p]={"data": [], "deviation": []}
                    data['indel'][p]={"data": [], "deviation": []}
                    
                for i in range(len(aux)):
                    mutsa={}
                    mutsb={}
                    for p in proteins:
                        mutsa[p]=[]
                        mutsb[p]=[]
                    
                    try:
                        if( len(str(aux.iloc[i,3]).split("-")) > 2 ):
                            #dat = datetime.fromisoformat(aux.iloc[i,3])
                            dat=aux.iloc[i,3]
                            if( int(dat.strftime("%Y")) > 2019):
                                
                                if(aux.iloc[i,14]!="NA" and aux.iloc[i,14]!=""):
                                    try:
                                        mutations=aux.iloc[i,14].replace("(","").replace(")","").split(",")
                                   
                                        for m in mutations:
                                        
                                            pr=m.split("_")[0]
                                            m=m.split("_")[1]
                                            if(m.find("del")!=-1 or m.find("ins")!=-1):
                                                mutsa[pr].append(m)
                                            else:
                                                mutsb[pr].append(m)
                                    except:
                                        pass
                                        
                                    for pr in proteins:
                                        indel[pr].append( len(mutsa[pr]) )
                                        mis[pr].append( len(mutsb[pr]) )
                                else:
                                    for pr in proteins:
                                        indel[pr].append(0)
                                        mis[pr].append(0)
            
                                per=dat.strftime(period)
            
                                if(ant!=per):
                                    if(ant!=""):
                                        labels.append(ant)
                                        for p in proteins:
                                            if( len(indel[p])==0 ):
                                                indel[p].append(0)
                                            if( len(indel[p])==1 ):
                                                indel[p].append(0)
                                                
                                            mean = statistics.median(indel[p])
                                            st = statistics.stdev(indel[p])
                                            data['indel'][p]["data"].append( mean )
                                            data['indel'][p]["deviation"].append( st )
                                            
                                            if( len(mis[p])==0 ):
                                                mis[p].append(0)
                                            if( len(mis[p])==1 ):
                                                mis[p].append(0)
                                                
                                            mean = statistics.median(mis[p])
                                            st = statistics.stdev(mis[p])
                                            data['missense'][p]["data"].append( mean )
                                            data['missense'][p]["deviation"].append( st )
                                    ant=per
                                    for p in proteins:
                                        indel[p]=[]
                                        mis[p]=[]
            
                    except:
                        pass
                    
                sum_=0
                for p in proteins:
                    sum_+=len(indel[p])
                
                if(sum_ > 0):
                    labels.append(ant)
                    for p in proteins:
                        if( len(indel[p])==0 ):
                            indel[p].append(0)
                        if( len(indel[p])==1 ):
                            indel[p].append(0)
                            
                        mean = statistics.median(indel[p])
                        st = statistics.stdev(indel[p])
                        data['indel'][p]["data"].append( mean )
                        data['indel'][p]["deviation"].append( st )
                        
                sum_=0
                for p in proteins:
                    sum_+=len(mis[p])
                
                if(sum_ > 0):
                    for p in proteins:
                        if( len(mis[p])==0 ):
                            mis[p].append(0)
                        if( len(mis[p])==1 ):
                            mis[p].append(0)
                            
                        mean = statistics.median(mis[p])
                        st = statistics.stdev(mis[p])
                        data['missense'][p]["data"].append( mean )
                        data['missense'][p]["deviation"].append( st )
            
                data['labels'] = labels 
                with open('data/mutation_indel/'+l.replace(" ","-")+"_"+names[co]+'_mutation_indel.json', 'w') as fp:
                    json.dump(data, fp)
            
            aux=df
            aux['Collection date'] = pd.to_datetime(aux["Collection date"])
            aux=aux.sort_values(by='Collection date')
            l="All"
            ant=""
            labels=[]
            
            indel={}
            mis={}
            for p in proteins:
                indel[p]=[]
                mis[p]=[]
                
                data['missense'][p]={"data": [], "deviation": []}
                data['indel'][p]={"data": [], "deviation": []}
                
            for i in range(len(aux)):
                mutsa={}
                mutsb={}
                for p in proteins:
                    mutsa[p]=[]
                    mutsb[p]=[]
                
                try:
                    if( len(str(aux.iloc[i,3]).split("-")) > 2 ):
                        #dat = datetime.fromisoformat(aux.iloc[i,3])
                        dat=aux.iloc[i,3]
                        if( int(dat.strftime("%Y")) > 2019):
                            
                            if(aux.iloc[i,14]!="NA" and aux.iloc[i,14]!=""):
                                try:
                                    mutations=aux.iloc[i,14].replace("(","").replace(")","").split(",")
                               
                                    for m in mutations:
                                    
                                        pr=m.split("_")[0]
                                        m=m.split("_")[1]
                                        if(m.find("del")!=-1 or m.find("ins")!=-1):
                                            mutsa[pr].append(m)
                                        else:
                                            mutsb[pr].append(m)
                                except:
                                    pass
                                    
                                for pr in proteins:
                                    indel[pr].append( len(mutsa[pr]) )
                                    mis[pr].append( len(mutsb[pr]) )
                            else:
                                for pr in proteins:
                                    indel[pr].append(0)
                                    mis[pr].append(0)
        
                            per=dat.strftime(period)
        
                            if(ant!=per):
                                if(ant!=""):
                                    labels.append(ant)
                                    for p in proteins:
                                        if( len(indel[p])==0 ):
                                            indel[p].append(0)
                                        if( len(indel[p])==1 ):
                                            indel[p].append(0)
                                            
                                        mean = statistics.median(indel[p])
                                        st = statistics.stdev(indel[p])
                                        data['indel'][p]["data"].append( mean )
                                        data['indel'][p]["deviation"].append( st )
                                        
                                        if( len(mis[p])==0 ):
                                            mis[p].append(0)
                                        if( len(mis[p])==1 ):
                                            mis[p].append(0)
                                            
                                        mean = statistics.median(mis[p])
                                        st = statistics.stdev(mis[p])
                                        data['missense'][p]["data"].append( mean )
                                        data['missense'][p]["deviation"].append( st )
                                ant=per
                                for p in proteins:
                                    indel[p]=[]
                                    mis[p]=[]
        
                except:
                    pass
                
            sum_=0
            for p in proteins:
                sum_+=len(indel[p])
            
            if(sum_ > 0):
                labels.append(ant)
                for p in proteins:
                    if( len(indel[p])==0 ):
                        indel[p].append(0)
                    if( len(indel[p])==1 ):
                        indel[p].append(0)
                        
                    mean = statistics.median(indel[p])
                    st = statistics.stdev(indel[p])
                    data['indel'][p]["data"].append( mean )
                    data['indel'][p]["deviation"].append( st )
                    
            sum_=0
            for p in proteins:
                sum_+=len(mis[p])
            
            if(sum_ > 0):
                for p in proteins:
                    if( len(mis[p])==0 ):
                        mis[p].append(0)
                    if( len(mis[p])==1 ):
                        mis[p].append(0)
                        
                    mean = statistics.median(mis[p])
                    st = statistics.stdev(mis[p])
                    data['missense'][p]["data"].append( mean )
                    data['missense'][p]["deviation"].append( st )
        
            data['labels'] = labels 
            with open("data/mutation_indel/All_"+names[co]+'_mutation_indel.json', 'w') as fp:
                json.dump(data, fp)
                
            co+=1
            
    def _run_forecasting(self, fil):
        df=pd.read_csv(fil, sep="\t")
        dfaux=pd.DataFrame()
        dfaux['ds']=df['date']
        dfaux['y']=df['y']
        gm_prophet = fbprophet.Prophet(changepoint_prior_scale=0.15)
        gm_prophet.fit(dfaux)
        # Make a future dataframe for 2 years
        gm_forecast = gm_prophet.make_future_dataframe(periods=30, freq='D')
        # Make predictions
        gm_forecast = gm_prophet.predict(gm_forecast)
        gm_prophet.plot(gm_forecast, xlabel = 'Date', ylabel = 'Mean Mutations')
        plt.title('Mutations - Predicting next days')
        plt.savefig(fil.replace('.tsv','')+".png")
        
    def get_forecasting_genome_data_general_location(self):
        locs=[]
        f=open("data/location.tsv","r")
        for line in f:
            locs.append(line.replace("\n",""))
        f.close()
        
        df=pd.read_csv("data/dates_days_mutations_all.txt", sep="\t")
        df.columns=['Collection date','Location', 'Muts']
        for i in range(len(df)):
             if( len(str(df.iloc[i,0])) == 8 ):
                 df.iloc[i,0] = str(df.iloc[i,0])[0:4]+"-"+str(df.iloc[i,0])[4:6]+"-"+str(df.iloc[i,0])[6:8]
             if( len(str(df.iloc[i,0])) == 7 ):
                 df.iloc[i,0] = str(df.iloc[i,0])[0:4]+"-"+str(df.iloc[i,0])[4:6]+"-0"+str(df.iloc[i,0])[6:8]
                 
        co=0
        #for period in ["w%V/%Y", "%m/%Y"]:
        period="%Y-%m-%d"
        
        for l in locs:
            f=open('data/forecasting/'+l.replace(" ","-")+'_forecasting_genome.tsv', 'w')
            f.write("date\ty\n")
            f.close()
            
            #if(not os.path.isfile('data/forecasting/'+l.replace(" ","-")+'_forecasting.tsv')):
            aux = df[ df['Location'].str.contains(l, na=False) ]
            aux['Collection date'] = pd.to_datetime(aux["Collection date"])
            aux=aux.sort_values(by='Collection date')
            ant=""
            means=[]
                
            for i in range(len(aux)):
                
                #try:
                if( len(str(aux.iloc[i,0])) > 7 ):
                    #dat = datetime.fromisoformat(aux.iloc[i,3])
                    dat=aux.iloc[i,0]
                    if( int(dat.strftime("%Y")) > 2019):
                        
                        m = aux.iloc[i,2]
                        means.append(m)
                        
                        per=dat.strftime(period)
    
                        if(ant!=per):
                            if(ant!=""):
                                mean = statistics.median(means)
                                with open('data/forecasting/'+l.replace(" ","-")+'_forecasting_genome.tsv', 'a') as gf:
                                    gf.write("%s\t%.2f\n" %(per, mean) )
                                    
                            means=[]
                            ant=per
    
                #except:
                #    pass
            
            if(len(means)>0):
                mean = statistics.median(means)
                with open('data/forecasting/'+l.replace(" ","-")+'_forecasting_genome.tsv', 'a') as gf:
                    gf.write("%s\t%.2f\n" %(per, mean) )
        
            #data['labels'] = labels 
            #with open('data/mutation_indel/'+l.replace(" ","-")+"_"+names[co]+'_mutation_indel.json', 'w') as fp:
            #    json.dump(data, fp)
        
        f=open('data/forecasting/All_forecasting_genome.tsv', 'w')
        f.write("date\ty\n")
        f.close()
        
        aux=df
        aux['Collection date'] = pd.to_datetime(aux["Collection date"])
        aux=aux.sort_values(by='Collection date')
        l="All"
        ant=""
        means=[]
            
        for i in range(len(aux)):
            
            #try:
            if( len(str(aux.iloc[i,0])) > 7 ):
                #dat = datetime.fromisoformat(aux.iloc[i,3])
                dat=aux.iloc[i,0]
                if( int(dat.strftime("%Y")) > 2019):
                    m = aux.iloc[i,2]
                    means.append(m)
                    
                    per=dat.strftime(period)

                    if(ant!=per):
                        if(ant!=""):
                            mean = statistics.median(means)
                            with open('data/forecasting/All_forecasting_genome.tsv', 'a') as gf:
                                gf.write("%s\t%.2f\n" %(per, mean) )
                                
                        means=[]
                        ant=per

            #except:
            #    pass
            
        if(len(means)>0):
                mean = statistics.median(means)
                with open('data/forecasting/All_forecasting_genome.tsv', 'a') as gf:
                    gf.write("%s\t%.2f\n" %(per, mean) )
        
    def get_forecasting_data_general_location(self):
        locs=[]
        f=open("data/location.tsv","r")
        for line in f:
            locs.append(line.replace("\n",""))
        f.close()
        
        proteins=[]
        f=open("data/proteins.tsv","r")
        for line in f:
            proteins.append(line.replace("\n",""))
        f.close()
        
        df=pd.read_csv("data/metadata.tsv", sep="\t")
        names=["weekly","monthly"]
        co=0
        #for period in ["w%V/%Y", "%m/%Y"]:
        period="%Y-%m-%d"
        data={"indel": {}, "missense": {}}
        
        for l in locs:
            f=open('data/forecasting/'+l.replace(" ","-")+'_forecasting.tsv', 'w')
            f.write("protein\tdate\ty\n")
            f.close()
            
            #if(not os.path.isfile('data/forecasting/'+l.replace(" ","-")+'_forecasting.tsv')):
            aux = df[ df['Location'].str.contains(l, na=False) ]
            aux['Collection date'] = pd.to_datetime(aux["Collection date"])
            aux=aux.sort_values(by='Collection date')
            ant=""
            labels=[]
            
            indel={}
            mis={}
            for p in proteins:
                indel[p]=[]
                mis[p]=[]
                
                data['missense'][p]={"data": [], "deviation": []}
                data['indel'][p]={"data": [], "deviation": []}
                
            for i in range(len(aux)):
                mutsa={}
                mutsb={}
                for p in proteins:
                    mutsa[p]=[]
                    mutsb[p]=[]
                
                #try:
                if( len(str(aux.iloc[i,3]).split("-")) > 2 ):
                    #dat = datetime.fromisoformat(aux.iloc[i,3])
                    dat=aux.iloc[i,3]
                    if( int(dat.strftime("%Y")) > 2019):
                        
                        if(aux.iloc[i,14]!="NA" and aux.iloc[i,14]!=""):
                            try:
                                mutations=aux.iloc[i,14].replace("(","").replace(")","").split(",")
                           
                                for m in mutations:
                                
                                    pr=m.split("_")[0]
                                    m=m.split("_")[1]
                                    if(m.find("del")!=-1 or m.find("ins")!=-1):
                                        mutsa[pr].append(m)
                                    else:
                                        mutsb[pr].append(m)
                            except:
                                pass
                                
                            for pr in proteins:
                                indel[pr].append( len(mutsa[pr]) )
                                mis[pr].append( len(mutsb[pr]) )
                        else:
                            for pr in proteins:
                                indel[pr].append(0)
                                mis[pr].append(0)
    
                        per=dat.strftime(period)
    
                        if(ant!=per):
                            if(ant!=""):
                                labels.append(ant)
                                for p in proteins:
                                    if( len(indel[p])==0 ):
                                        indel[p].append(0)
                                    if( len(indel[p])==1 ):
                                        indel[p].append(0)
                                        
                                    mean = statistics.median(indel[p])
                                    st = statistics.stdev(indel[p])
                                    data['indel'][p]["data"].append( mean )
                                    data['indel'][p]["deviation"].append( st )
                                    
                                    if( len(mis[p])==0 ):
                                        mis[p].append(0)
                                    if( len(mis[p])==1 ):
                                        mis[p].append(0)
                                        
                                    mean = statistics.median(mis[p])
                                    st = statistics.stdev(mis[p])
                                    data['missense'][p]["data"].append( mean )
                                    data['missense'][p]["deviation"].append( st )
                                    
                                    with open('data/forecasting/'+l.replace(" ","-")+'_forecasting.tsv', 'a') as gf:
                                        gf.write("%s\t%s\t%.2f\n" %(p, per, mean) )
                                    
                            ant=per
                            for p in proteins:
                                indel[p]=[]
                                mis[p]=[]
    
                #except:
                #    pass
                
            sum_=0
            for p in proteins:
                sum_+=len(indel[p])
            
            if(sum_ > 0):
                labels.append(ant)
                for p in proteins:
                    if( len(indel[p])==0 ):
                        indel[p].append(0)
                    if( len(indel[p])==1 ):
                        indel[p].append(0)
                        
                    mean = statistics.median(indel[p])
                    st = statistics.stdev(indel[p])
                    data['indel'][p]["data"].append( mean )
                    data['indel'][p]["deviation"].append( st )
                    
            sum_=0
            for p in proteins:
                sum_+=len(mis[p])
            
            if(sum_ > 0):
                for p in proteins:
                    if( len(mis[p])==0 ):
                        mis[p].append(0)
                    if( len(mis[p])==1 ):
                        mis[p].append(0)
                        
                    mean = statistics.median(mis[p])
                    st = statistics.stdev(mis[p])
                    data['missense'][p]["data"].append( mean )
                    data['missense'][p]["deviation"].append( st )
                                    
                    with open('data/forecasting/'+l.replace(" ","-")+'_forecasting.tsv', 'a') as gf:
                        gf.write("%s\t%s\t%.2f\n" %(p, per, mean) )
        
            #data['labels'] = labels 
            #with open('data/mutation_indel/'+l.replace(" ","-")+"_"+names[co]+'_mutation_indel.json', 'w') as fp:
            #    json.dump(data, fp)
        
        f=open('data/forecasting/forecasting.tsv', 'w')
        f.write("protein\tdate\ty\n")
        f.close()
        
        aux=df
        aux['Collection date'] = pd.to_datetime(aux["Collection date"])
        aux=aux.sort_values(by='Collection date')
        l="All"
        ant=""
        labels=[]
        
        indel={}
        mis={}
        for p in proteins:
            indel[p]=[]
            mis[p]=[]
            
            data['missense'][p]={"data": [], "deviation": []}
            data['indel'][p]={"data": [], "deviation": []}
            
        for i in range(len(aux)):
            mutsa={}
            mutsb={}
            for p in proteins:
                mutsa[p]=[]
                mutsb[p]=[]
            
            #try:
            if( len(str(aux.iloc[i,3]).split("-")) > 2 ):
                #dat = datetime.fromisoformat(aux.iloc[i,3])
                dat=aux.iloc[i,3]
                if( int(dat.strftime("%Y")) > 2019):
                    
                    if(aux.iloc[i,14]!="NA" and aux.iloc[i,14]!=""):
                        try:
                            mutations=aux.iloc[i,14].replace("(","").replace(")","").split(",")
                       
                            for m in mutations:
                            
                                pr=m.split("_")[0]
                                m=m.split("_")[1]
                                if(m.find("del")!=-1 or m.find("ins")!=-1):
                                    mutsa[pr].append(m)
                                else:
                                    mutsb[pr].append(m)
                        except:
                            pass
                            
                        for pr in proteins:
                            indel[pr].append( len(mutsa[pr]) )
                            mis[pr].append( len(mutsb[pr]) )
                    else:
                        for pr in proteins:
                            indel[pr].append(0)
                            mis[pr].append(0)

                    per=dat.strftime(period)

                    if(ant!=per):
                        if(ant!=""):
                            labels.append(ant)
                            for p in proteins:
                                if( len(indel[p])==0 ):
                                    indel[p].append(0)
                                if( len(indel[p])==1 ):
                                    indel[p].append(0)
                                    
                                mean = statistics.median(indel[p])
                                st = statistics.stdev(indel[p])
                                data['indel'][p]["data"].append( mean )
                                data['indel'][p]["deviation"].append( st )
                                
                                if( len(mis[p])==0 ):
                                    mis[p].append(0)
                                if( len(mis[p])==1 ):
                                    mis[p].append(0)
                                    
                                mean = statistics.median(mis[p])
                                st = statistics.stdev(mis[p])
                                data['missense'][p]["data"].append( mean )
                                data['missense'][p]["deviation"].append( st )
                                    
                                with open('data/forecasting/All_forecasting.tsv', 'a') as gf:
                                    gf.write("%s\t%s\t%.2f\n" %(p, per, mean) )
                        ant=per
                        for p in proteins:
                            indel[p]=[]
                            mis[p]=[]

            #except:
            #    pass
            
        sum_=0
        for p in proteins:
            sum_+=len(indel[p])
        
        if(sum_ > 0):
            labels.append(ant)
            for p in proteins:
                if( len(indel[p])==0 ):
                    indel[p].append(0)
                if( len(indel[p])==1 ):
                    indel[p].append(0)
                    
                mean = statistics.median(indel[p])
                st = statistics.stdev(indel[p])
                data['indel'][p]["data"].append( mean )
                data['indel'][p]["deviation"].append( st )
                
        sum_=0
        for p in proteins:
            sum_+=len(mis[p])
        
        if(sum_ > 0):
            for p in proteins:
                if( len(mis[p])==0 ):
                    mis[p].append(0)
                if( len(mis[p])==1 ):
                    mis[p].append(0)
                    
                mean = statistics.median(mis[p])
                st = statistics.stdev(mis[p])
                data['missense'][p]["data"].append( mean )
                data['missense'][p]["deviation"].append( st )
                                    
                with open('data/forecasting/All_forecasting.tsv', 'a') as gf:
                    gf.write("%s\t%s\t%.2f\n" %(p, per, mean) )
                    
        for f in os.listdir('data/forecasting/'):
            if(f.endswith(".tsv")):
                self._run_forecasting('data/forecasting/'+f)
            
    def get_conditional_probabilities(self):
        groups={ 'Non Polar': ['G','A','V','C','P','L', 'I', 'M','W','F'], 'Polar': ['S','T','Y','N','Q'], 'Positive Charge': ['K','R','H'], 'Negative Charge': ['D','E'] }
    
        f=open("data/report_position_probabilities_by_lineage.tsv","w")
        f.write("Lineage\tProtein\tPosition\tEntropy\tAA Group Change\tAA Change\tCount\tProportion\n")
        f.close()
        
        dd=pd.read_csv("data/metadata.tsv", sep="\t")
        vocs=["P.1","B.1.1.7","B.1.351","B.1.617.2"]
        for m in vocs:
            probs={}
            aux=list(dd[ dd['Pango lineage']==m ]['AA Substitutions'].values)
            for c in aux:
                try:
                    lin="NA"
                    if(c!="" and c!="unknown" and c!="nan"):
                        lin=c.replace("(","").replace(")","").split(",")
                        
                    if(lin!="NA" and lin!="None"):
                        for l in lin:
                            protein=l.split("_")[0]
                            mut=l.split("_")[1]
                            
                            if(not mut.startswith("ins")):
                                if(not protein in probs.keys()):
                                    probs[protein]={}
                                    
                                end=-1
                                ref=mut[0]
                                change=mut[end]
                                
                                if(mut.endswith("del")):
                                    end=-3
                                    change="del"
                                if(mut.endswith("stop")):
                                    end=-4
                                    change="stop"
                                
                                pos=mut[1:end]
                                try:
                                    ps=int(pos)
                                    if(not pos in probs[protein].keys()):
                                        probs[protein][pos]={}
                                        
                                    if(not ref+"-"+change in probs[protein][pos].keys()):
                                        probs[protein][pos][ref+"-"+change]=0
                                    probs[protein][pos][ref+"-"+change]+=1
                                except:
                                    pass
                except:
                    pass
            
            for pr in probs.keys():
                for p in probs[pr].keys():
                    entropy=0
                    
                    nabsolute=[]
                    nrelative=[]   
                    aas=probs[pr][p].keys()
                    total=sum(probs[pr][p].values())
                    
                    group_changes=[]
                    for aa in probs[pr][p].keys():
                        ref=aa.split("-")[0]
                        change=aa.split("-")[1]
                        refgr=""
                        for gr in groups.keys():
                            if(ref in groups[gr]):
                                refgr=gr
                                break
                        
                        chagr=""
                        for gr in groups.keys():
                            if(change in groups[gr]):
                                chagr=gr
                                break
                        group_changes.append(refgr+" to "+chagr)
                        
                        count=probs[pr][p][aa]
                        paa=count/total
                        
                        nabsolute.append( str(count) )
                        nrelative.append( str(paa) )
                        
                        entropy+=paa*(math.log(paa, 2))
                    
                    entropy *= -1
                    
                    c=0
                    for a in aas:
                        with open("data/report_position_probabilities_by_lineage.tsv","a") as gf:
                            gf.write("%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\n" %(m, pr, p, entropy, group_changes[c], a.replace("-"," to "), nabsolute[c], nrelative[c] ) )
                        c+=1
    
    def get_conditional_probabilities_by_month(self):
        locs=[]
        f=open("data/location.tsv","r")
        for line in f:
            locs.append(line.replace("\n",""))
        f.close()
        
        df=pd.read_csv("data/metadata.tsv", sep="\t")
        mutations=df['AA Substitutions'].values
        locations=df['Location'].values
        c=0
        for l in locations:
            cot=""
            for lc in locs:
                if (l.find(lc)!=-1):
                    cot=lc
                    break
                
            locations[c]=cot
            c+=1
            
        groups={ 'Non Polar': ['G','A','V','C','P','L', 'I', 'M','W','F'], 'Polar': ['S','T','Y','N','Q'], 'Positive Charge': ['K','R','H'], 'Negative Charge': ['D','E'] }
    
        f=open("data/report_position_probabilities_by_lineage_period.tsv","w")
        f.write("Lineage\tMonth\tYear\tProtein\tPosition\tEntropy\tAA Group Change\tAA Change\tCount\tProportion\n")
        f.close()
        
        dd=pd.read_csv("data/metadata.tsv", sep="\t")
        vocs=["P.1","B.1.1.7","B.1.351","B.1.617.2"]
        for m in vocs:
            probs={}
            
            aux=dd[ dd['Pango lineage']==m ]
            #dates = pd.to_datetime(aux["Collection date"])
            #aux["Collection date"]=dates
            #aux=aux.sort("Collection date")
            #aux=aux.sort_values(by=['Collection date'])
            mutations=list(aux['AA Substitutions'].values)
            collection=list(aux['Collection date'].values)
            
            period="%m/%Y"
            
            count=-1
            ant=""
            for c in mutations:
                count+=1
                try:
                    dat=collection[count]
                    
                    dat = date.fromisoformat(dat)
                    per=dat.strftime(period)
                    print(per)
                    if( int(dat.strftime("%Y")) > 2019):
                        lin="NA"
                        if(c!="" and c!="unknown" and c!="nan"):
                            lin=c.replace("(","").replace(")","").split(",")
                            
                        if(lin!="NA" and lin!="None"):
                            if(not per in probs.keys()):
                                probs[per]={}
                                
                            for l in lin:
                                protein=l.split("_")[0]
                                mut=l.split("_")[1]
                                
                                if(not mut.startswith("ins")):
                                    if(not protein in probs[per].keys()):
                                        probs[per][protein]={}
                                        
                                    end=-1
                                    change=mut[end]
                                    ref=mut[0]
                                        
                                    if(mut.endswith("del")):
                                        end=-3
                                        change="del"
                                    if(mut.endswith("stop")):
                                        end=-4
                                        change="stop"
                                    
                                    pos=mut[1:end]
                                    try:
                                        ps=int(pos)
                                        if(not pos in probs[per][protein].keys()):
                                            probs[per][protein][pos]={}
                                            
                                        if(not ref+"-"+change in probs[per][protein][pos].keys()):
                                            probs[per][protein][pos][ref+"-"+change]=0
                                        probs[per][protein][pos][ref+"-"+change]+=1
                                    except:
                                        pass
                except:
                    pass
            
            for per in probs.keys():
                for pr in probs[per].keys():
                    for p in probs[per][pr].keys():
                        entropy=0
                        
                        nabsolute=[]
                        nrelative=[]   
                        aas=probs[per][pr][p].keys()
                        total=sum(probs[per][pr][p].values())
                        
                        group_changes=[]
                        for aa in probs[per][pr][p].keys():
                            ref=aa.split("-")[0]
                            change=aa.split("-")[1]
                            refgr=""
                            for gr in groups.keys():
                                if(ref in groups[gr]):
                                    refgr=gr
                                    break
                            
                            chagr=""
                            for gr in groups.keys():
                                if(change in groups[gr]):
                                    chagr=gr
                                    break
                            group_changes.append(refgr+" to "+chagr)
                                
                            count=probs[per][pr][p][aa]
                            paa=count/total
                            
                            nabsolute.append( str(count) )
                            nrelative.append( str(paa) )
                            
                            entropy+=paa*(math.log(paa, 2))
                        
                        entropy *= -1
                    
                        c=0
                        for a in aas:
                            month=per.split("/")[0]
                            if(month.startswith("0")):
                                month=month[1:]
                                
                            year=per.split("/")[1]
                            
                            with open("data/report_position_probabilities_by_lineage_period.tsv","a") as gf:
                                gf.write("%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\n" %(m, month, year, pr, p, entropy, group_changes[c], a.replace("-"," to "), nabsolute[c], nrelative[c] ) )
                            
                            c+=1
            
        for loc in locs:
            f=open("data/entropy/"+loc.replace(" ","-")+"_report_position_probabilities_by_lineage_period.tsv","w")
            f.write("Lineage\tMonth\tYear\tProtein\tPosition\tEntropy\tAA Group Change\tAA Change\tCount\tProportion\n")
            f.close()
            
            for m in vocs:
                probs={}
            
                aux=dd[ dd['Pango lineage']==m  ]
                aux=aux[ aux['Location'].str.contains(loc, na=False) ]
                #dates = pd.to_datetime(aux["Collection date"])
                #aux["Collection date"]=dates
                #aux=aux.sort("Collection date")
                #aux=aux.sort_values(by=['Collection date'])
                mutations=list(aux['AA Substitutions'].values)
                collection=list(aux['Collection date'].values)
                
                period="%m/%Y"
                
                count=-1
                ant=""
                for c in mutations:
                    count+=1
                    try:
                        dat=collection[count]
                        
                        dat = date.fromisoformat(dat)
                        per=dat.strftime(period)
                        print(per)
                        if( int(dat.strftime("%Y")) > 2019):
                            lin="NA"
                            if(c!="" and c!="unknown" and c!="nan"):
                                lin=c.replace("(","").replace(")","").split(",")
                                
                            if(lin!="NA" and lin!="None"):
                                if(not per in probs.keys()):
                                    probs[per]={}
                                    
                                for l in lin:
                                    protein=l.split("_")[0]
                                    mut=l.split("_")[1]
                                    
                                    if(not mut.startswith("ins")):
                                        if(not protein in probs[per].keys()):
                                            probs[per][protein]={}
                                            
                                        end=-1
                                        change=mut[end]
                                        ref=mut[0]
                                            
                                        if(mut.endswith("del")):
                                            end=-3
                                            change="del"
                                        if(mut.endswith("stop")):
                                            end=-4
                                            change="stop"
                                        
                                        pos=mut[1:end]
                                        try:
                                            ps=int(pos)
                                            if(not pos in probs[per][protein].keys()):
                                                probs[per][protein][pos]={}
                                                
                                            if(not ref+"-"+change in probs[per][protein][pos].keys()):
                                                probs[per][protein][pos][ref+"-"+change]=0
                                            probs[per][protein][pos][ref+"-"+change]+=1
                                        except:
                                            pass
                    except:
                        pass
                
                for per in probs.keys():
                    for pr in probs[per].keys():
                        for p in probs[per][pr].keys():
                            entropy=0
                            
                            nabsolute=[]
                            nrelative=[]   
                            aas=probs[per][pr][p].keys()
                            total=sum(probs[per][pr][p].values())
                            
                            group_changes=[]
                            for aa in probs[per][pr][p].keys():
                                ref=aa.split("-")[0]
                                change=aa.split("-")[1]
                                refgr=""
                                for gr in groups.keys():
                                    if(ref in groups[gr]):
                                        refgr=gr
                                        break
                                
                                chagr=""
                                for gr in groups.keys():
                                    if(change in groups[gr]):
                                        chagr=gr
                                        break
                                group_changes.append(refgr+" to "+chagr)
                                    
                                count=probs[per][pr][p][aa]
                                paa=count/total
                                
                                nabsolute.append( str(count) )
                                nrelative.append( str(paa) )
                                
                                entropy+=paa*(math.log(paa, 2))
                            
                            entropy *= -1
                        
                            c=0
                            for a in aas:
                                month=per.split("/")[0]
                                if(month.startswith("0")):
                                    month=month[1:]
                                    
                                year=per.split("/")[1]
                                    
                                with open("data/entropy/"+loc.replace(" ","-")+"_report_position_probabilities_by_lineage_period.tsv","a") as gf:
                                    gf.write("%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\n" %(m, month, year, pr, p, entropy, group_changes[c], a.replace("-"," to "), nabsolute[c], nrelative[c] ) )
                                
                                c+=1
                        
    def process_info_domain(self):
        domain_info={}
        f=open("data/prosite.dat","r")
        for line in f:
            l=line.replace('\n','')
            if(l.startswith("ID")):
                id_=l.split("   ")[1].split(";")[0]
                
            if(l.startswith("AC")):
                ac=l.split("   ")[1].split(";")[0]
                domain_info[id_]={'prosite': ac, 'description': ""}
                
            if(l.startswith("DE")):
                de=l.split("   ")[1].split(";")[0]
                domain_info[id_]['description']=de
        f.close()
        
        with open('data/info_domains.json', 'w') as fp:
            json.dump(domain_info, fp)
        
    def analyze_domains(self):
        domains={}
        # use included_sequences_2
        
        f=open("../../pipeline_mutation_vigilance/included_sequences_2.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            
            protein=l[1]
            
            mutaions_list=l[3].split(";")
            
            domains_list=l[5].split(";")
            coords=l[6].split(";")
            
            
            """for m in mutaions_list:
                if(m!="NA"):
                    pos=int(m[1:-1])
                    a=0
                    for c in coords:
                        if(not domains_list[a] in mutsd):
                            mutsd[domains_list[a]]=[]
                        
                        if(c!="NA"):
                            start=int(c.split("-")[0])
                            end=int(c.split("-")[1])
                            if(pos>=start and pos<=end):
                                if(not m in mutsd[domains_list[a]]):
                                    mutsd[domains_list[a]].append(m)
                        a+=1"""
                        
            if(not protein in domains.keys()):
                domains[protein]={}
            
            i=0
            for d in domains_list:
                if(d!="NA"):
                    if(not d in domains[protein]):
                        domains[protein][d]=[0, 0, []]
                    domains[protein][d][0]+=1
                    
                    c=coords[i]
                    for m in mutaions_list:
                        if(m!="NA"):
                            pos=int(m[1:-1])
                            start=int(c.split("-")[0])
                            end=int(c.split("-")[1])
                            if(pos>=start and pos<=end):
                                if(not m in domains[protein][d][2]):
                                    domains[protein][d][2].append(m)
                                    domains[protein][d][1]+=1
                i+=1
        f.close()
        
        # table: protein, domain, occurrences, mutations, mutation_occurrences
        
        f=open("data/report_domains.tsv","w")
        f.write("label\tx\ty\tqtd_mutations\tmutations\n")
        for pr in domains.keys():
            print(domains[pr].keys())
            for d in domains[pr].keys():
                f.write("%s\t%s\t%i\t%i\t%s\n" %( pr, d, domains[pr][d][0], domains[pr][d][1], ';'.join(domains[pr][d][2]) ) )
        f.close()
    
    def lineageCount_by_mutation(self):
        df=pd.read_csv("data/metadata.tsv", sep="\t")
        mutations=df['AA Substitutions'].values
        lineage=df['Pango lineage'].values
        c=0
        counts={}
        for muts in mutations:
            try:
                lin=lineage[c]
                
                if(muts!='' and muts!='NA'):
                    muts=muts.replace("(","").replace(")","").split(",")
                    for m in muts:
                        pr=m.split("_")[0]
                        mut=m.split("_")[1]
                        if(not m in counts.keys()):
                            counts[m]={}
                        
                        if(lin!='' and lin!='unknown'):
                            if(not lin in counts[m].keys()):
                                counts[m][lin]=0
                            
                            counts[m][lin]+=1
            except:
                pass
            
            c+=1
            
        with open('data/lineageCount_by_mutation.json', 'w') as fp:
                json.dump(counts, fp)
                
    def countryCount_by_mutation(self):
        locs=[]
        f=open("data/location.tsv","r")
        for line in f:
            locs.append(line.replace("\n",""))
        f.close()
        
        df=pd.read_csv("data/metadata.tsv", sep="\t")
        mutations=df['AA Substitutions'].values
        locations=df['Location'].values
        c=0
        for l in locations:
            cot=""
            for lc in locs:
                if (l.find(lc)!=-1):
                    cot=lc
                    break
                
            locations[c]=cot
            c+=1
            
        c=0
        counts={}
        for muts in mutations:
            try:
                lin=locations[c]
                
                if(muts!='' and muts!='NA'):
                    muts=muts.replace("(","").replace(")","").split(",")
                    for m in muts:
                        pr=m.split("_")[0]
                        mut=m.split("_")[1]
                        if(not m in counts.keys()):
                            counts[m]={}
                        
                        if(lin!='' and lin!='unknown'):
                            if(not lin in counts[m].keys()):
                                counts[m][lin]=0
                            
                            counts[m][lin]+=1
            except:
                pass
            
            c+=1
            
        with open('data/countryCount_by_mutation.json', 'w') as fp:
                json.dump(counts, fp)
    
    def mutationCount_by_lineage(self):
        
        df=pd.read_csv("data/metadata.tsv", sep="\t")
        mutations=df['AA Substitutions'].values
        locations=df['Pango lineage'].values
        
        c=0
        counts={}
        for muts in mutations:
            try:
                lin=locations[c]
                
                if(muts!='' and muts!='NA'):
                    muts=muts.replace("(","").replace(")","").split(",")
                    for m in muts:
                        pr=m.split("_")[0]
                        mut=m.split("_")[1]
                        
                        if(lin!='' and lin!='unknown'):
                            if(not lin in counts.keys()):
                                counts[lin]={}
                            
                            if(not m in counts[lin].keys()):
                                counts[lin][m]=0
                            counts[lin][m]+=1
                        
            except:
                pass
            
            c+=1
            
        with open('data/mutationCount_by_lineage.json', 'w') as fp:
                json.dump(counts, fp)
    
    def uniqueMutationCount_by_lineage(self):
        with open('data/mutationCount_by_lineage.json', 'r') as fp:
            dat1 = json.load(fp)
            
        with open('data/lineageCount_by_mutation.json', 'r') as fp:
            dat2 = json.load(fp)
        
        unique={}
        for m in dat2:
            lineages = list(dat2[m].keys())
            if( len(lineages) == 1):
                if(not lineages[0] in unique.keys()):
                    unique[lineages[0]] = {}
                    
                if(not m in unique[lineages[0]].keys()):
                    unique[lineages[0]][m] = 0
                    
                unique[lineages[0]][m] = dat1[lineages[0]][m]
            
        with open('data/uniqueMutationCount_by_lineage.json', 'w') as fp:
                json.dump(unique, fp)
    
    def mutationCount(self):
        df=pd.read_csv("data/metadata.tsv", sep="\t")
        mutations=df['AA Substitutions'].values
        
        c=0
        counts={}
        for muts in mutations:
            try:
                if(muts!='' and muts!='NA'):
                    muts=muts.replace("(","").replace(")","").split(",")
                    for m in muts:
                        pr=m.split("_")[0]
                        mut=m.split("_")[1]
                        if(not m in counts.keys()):
                            counts[m]=0
                            
                        counts[m]+=1
            except:
                pass
            
            c+=1
            
        with open('data/mutationCount.json', 'w') as fp:
                json.dump(counts, fp)
    
    def lineageCount(self):
        df=pd.read_csv("data/metadata.tsv", sep="\t")
        lins=df['Pango lineage'].values
        
        c=0
        counts={}
        for m in lins:
            try:
                if(not m in counts.keys()):
                    counts[m]=0
                    
                counts[m]+=1
            except:
                pass
            
            c+=1
            
        with open('data/lineageCount.json', 'w') as fp:
                json.dump(counts, fp)
                
    def locationsCount(self):
        locs=[]
        f=open("data/location.tsv","r")
        for line in f:
            locs.append(line.replace("\n",""))
        f.close()
        
        df=pd.read_csv("data/metadata.tsv", sep="\t")
        lins=df['Location'].values
        
        c=0
        counts={}
        for m in lins:
            try:
                cot=""
                for lc in locs:
                    if (m.find(lc)!=-1):
                        cot=lc
                        break
                    
                if(not cot in counts.keys()):
                    counts[cot]=0
                    
                counts[cot]+=1
            except:
                pass
            
            c+=1
            
        with open('data/locationsCount.json', 'w') as fp:
                json.dump(counts, fp)
                
    def run(self):
        
        #self.get_mutation_indel_data_general_location()
        #self.get_mutation_indel_data_lineage()
        #self.analyze_domains()
        #self.process_info_domain()
        
        #self.lineageCount_by_mutation()
        #self.countryCount_by_mutation()
        #self.mutationCount_by_lineage()
        #self.uniqueMutationCount_by_lineage()
        
        #self.get_conditional_probabilities()
        #self.get_conditional_probabilities_by_month()
        
        self.get_demography_data()
        
        #self.mutationCount()
        #self.lineageCount()
        #self.locationsCount()
        
        #self.get_all_snap_results()
        
        #self.get_forecasting_data_general_location()
        #self.get_forecasting_genome_data_general_location()
        
a=Generate_processed_data()
a.run()
