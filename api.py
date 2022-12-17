from flask import Flask
from flask import request, jsonify, send_from_directory, send_file
from flask_cors import CORS, cross_origin

app = Flask(__name__)
cors = CORS(app)

import pandas as pd
from datetime import datetime
import json
import os
import statistics
from faker import Factory

import fbprophet
from matplotlib import pyplot as plt

import simplejson as json

import uuid
import shutil
from auxiliar_functions import _get_mutations, _get_domain, _get_muts_in_domains, send_success_email

import Levenshtein

#Acesso a libs do C/C++
import ctypes

#directory="/var/www/html/"
#host="localhost:5000"
#ip="127.0.0.1"

directory="/var/www/html/foca_backend/"
host="127.0.0.1"
ip="127.0.0.1"

@app.route("/")
def principal():
    return "FOCA"

@app.route('/static/results/<path:filename>')
def serve_static(filename):
    root_dir = os.path.dirname(os.getcwd())
    root_dir = directory
    return send_from_directory( os.path.join(root_dir, 'data_exportation'), filename)

@app.route('/get_proteins', methods=['GET'])
def section0_get_proteins():
    proteins=[]
    f=open(directory+"data/proteins.tsv","r")
    for line in f:
        l=line.replace("\n", "")
        proteins.append(l)
    f.close()

    resp={"msg": proteins}

    return json.dumps(resp)

@app.route('/get_locations', methods=['GET'])
def section0_get_locations():
    locations=[]
    f=open(directory+"data/location.tsv","r")
    for line in f:
        l=line.replace("\n", "")
        try:
            a=int(l)
        except:
            wrong=['bat','canine','cat','dog','dog ','env','gorilla','hamster','leopard','lion','mink','monkey','mouse','pangolin','snow_leopard','tiger']
            if(not l in wrong):
                locations.append(l)
    f.close()

    resp={"msg": locations}

    return json.dumps(resp)

@app.route('/get_lineages', methods=['GET'])
def section0_get_lineages():
    lineages=[]
    f=open(directory+"data/lineage.tsv","r")
    for line in f:
        l=line.replace("\n", "")
        lineages.append(l)
    f.close()

    resp={"msg": lineages}

    return json.dumps(resp)

@app.route('/get_status', methods=['GET'])
def section0_get_status():
    status=[]
    f=open(directory+"data/demography/list_status.tsv","r")
    for line in f:
        l=line.replace("\n", "")
        status.append(l)
    f.close()

    resp={"msg": status}

    return json.dumps(resp)

@app.route('/last_update', methods=['GET'])
def section1_last_update():
    lu=""
    f=open(directory+"data/last_update.tsv","r")
    for line in f:
        l=line.replace("\n", "")
        lu=l
    f.close()

    resp={"msg": lu}
    return json.dumps(resp)

# Page demography analysis
@app.route('/get_demography_plots', methods=['GET'])
def section_get_plots_demography():
    resp={"error": ""}
    with open(directory+'data/demography/global_data.json', 'r') as fp:
        dat = json.load(fp)  
    resp["global_geo_data"]=dat
    
    with open(directory+'data/demography/global_info.json', 'r') as fp:
        dat = json.load(fp)  
    resp["global_geo_info"]=dat
    
    with open(directory+'data/demography/by_age_data.json', 'r') as fp:
        dat = json.load(fp)  
    resp["plot_age"]=dat
    
    with open(directory+'data/demography/by_voc_data.json', 'r') as fp:
        dat = json.load(fp)  
    resp["plot_voc"]=dat
    
    with open(directory+'data/demography/by_patient_status_data.json', 'r') as fp:
        dat = json.load(fp)  
    resp["plot_status"]=dat
    
    with open(directory+'data/demography/age_by_voc_data.json', 'r') as fp:
        dat = json.load(fp)  
    resp["plot_age_voc"]=dat
    
    with open(directory+'data/demography/age_by_status_data.json', 'r') as fp:
        dat = json.load(fp)  
    resp["plot_age_status"]=dat
    
    with open(directory+'data/demography/status_by_voc_data.json', 'r') as fp:
        dat = json.load(fp)  
    resp["plot_status_voc"]=dat

    return json.dumps(resp)

@app.route('/get_table_demography/<lineage>/<location>/<status>/<gender>/<age>', methods=['GET'])
def section_get_table_demography(lineage, location, status, gender, age):
    data={"error": "", "columns": [], "table": [], "file_export": "/static/results/"}
    
    uui=str(uuid.uuid4())
    
    df=pd.read_csv(directory+'data/demography/metadata_demography.tsv', sep='\t')
    conditions=[]
    lineage=lineage.replace("-",".")
    if(lineage!="All"):
        conditions.append(' ( df["Lineage"]=="'+lineage+'" ) ')
    
    if(location!="All"):
        conditions.append(' ( df["Location"]=="'+location+'" ) ')
    
    if(status!="All"):
        conditions.append(' ( df["Status"].str.contains("'+status+'", na=False) ) ')
    
    if(gender!="All"):
        conditions.append(' ( df["Gender"]=="'+gender+'" ) ')
    
    if(age!=""):
        conditions.append(' ( df["Age"]>'+age+' ) ')
        
    condition="( "+' & '.join(conditions)+" )"
    if(len(conditions)>0):
        df = df[ eval(condition) ]
        
    if(len(df)!=0):
        #g=open(directory+'data_exportation/'+uui+".tsv", "w")
        columns=[]
        data['columns']=[]
        for col in df.columns:
            columns.append(col)
            data['columns'].append({ "title" : col } )
        
        data['columns']=data['columns'][:-1]
        
        df.to_csv(directory+'data_exportation/'+uui+".tsv", sep="\t")
        df=df.iloc[:1000, :-1]
        #g.write( ('\t'.join(columns))+"\n")
        
        for i in range(len(df)):
            dat=list(df.iloc[i, :])
            c=0
            for d in dat:
                dat[c]=str(d)
                c+=1 
                
            if(i<1000):
                data["table"].append( dat )
            #g.write( ('\t'.join(dat))+"\n")
        #g.close()
        
        data['file_export'] += uui+'.tsv'
    else:
        data['error']='There are no matches for your filters.'
 
    return json.dumps(data)

# Page entropy vocs analysis
@app.route('/entropy_position_analysis/<lineage>/<protein>/<position>/<effect>/<count>/<proportion>', methods=['GET'])
def section_entropy_position_analysis(lineage, protein, position, effect, count, proportion):
    groups={ 'Non Polar': ['G','A','V','C','P','L', 'I', 'M','W','F'], 'Polar': ['S','T','Y','N','Q'], 'Positive Charge': ['K','R','H'], 'Negative Charge': ['D','E'] }
    revg={"del": "-"}
    for g in groups:
        for aa in groups[g]:
            revg[aa]=g
            
    uui=str(uuid.uuid4())
    
    sr={}
    with open(directory+'data/snap_data.json', 'r') as fp:
        sr = json.load(fp)  
        
    data={"error": "","table_entropy": [], "file_export": "/static/results/"}
    
    # Generating table
    df=pd.read_csv(directory+'data/report_position_probabilities_by_lineage.tsv', sep='\t')
    conditions=[]
    lineage=lineage.replace("-",".")
    proportion=proportion.replace("-",".")
    if(lineage!="All"):
        conditions.append(' ( df["Lineage"]=="'+lineage+'" ) ')
    
    if(protein!="All"):
        cproteins=[]
        for m in protein.split(','):
            cproteins.append(' ( df["Protein"]=="'+m+'" ) ')
        conditions.append("( "+' | '.join(cproteins)+" )")
        
    if(position!="All"):
        conditions.append(' ( df["Position"]=="'+position+'" ) ')
    
    if(count!=""):
        conditions.append(' ( df["Count"]>'+count+' ) ')
    
    if(proportion!=""):
        conditions.append(' ( df["Proportion"]>'+proportion+' ) ')
        
    condition="( "+' & '.join(conditions)+" )"
    if(len(conditions)>0):
        df = df[ eval(condition) ]
        
    if(len(df)!=0):
        g=open(directory+'data_exportation/'+uui+".tsv", "w")
        columns=[]
        data['columns']=[]
        for c in df.columns:
            columns.append(c)
            data['columns'].append({ "title" : c } )
        
        columns+=["SNAP effect","SNAP prediction"]
        data['columns']+=[{ "title" : 'SNAP effect' }, { "title": 'SNAP prediction' }]
        
        g.write( ('\t'.join(columns))+"\n")
        
        for i in range(len(df)):
            protein=df.iloc[i,1].lower()
            position = str(df.iloc[i,2])
            aas = df.iloc[i,5].split(" to ")
            
            sn=["-","-"]
            if( protein+"_"+aas[0]+position+aas[1] in sr.keys()):
                sn=sr[protein+"_"+aas[0]+position+aas[1]]
            
            if(sn[0]==effect or effect=="All"):
                dat=list(df.iloc[i, :])+sn
                if(i<1000):
                    data["table_entropy"].append( dat )
                
                c=0
                for d in dat:
                    dat[c]=str(d)
                    c+=1 
                g.write( ('\t'.join(dat))+"\n")
        g.close()
        
        data['file_export'] += uui+'.tsv'
    else:
        data['error']='There are no matches for your filters.'
    
    return json.dumps(data)

@app.route('/plot_entropy_position_analysis_by_period/<lineage>/<location>', methods=['GET'])
def section_plot_entropy_position_analysis_by_period(lineage, location):
    # Generating plot over time
    fake = Factory.create()
    
    resp={}
    
    if(location=="All"):
        df=pd.read_csv(directory+'data/report_position_probabilities_by_lineage_period.tsv', sep='\t')
    else:
        df=pd.read_csv(directory+"data/entropy/"+location.replace(" ", "-")+"_report_position_probabilities_by_lineage_period.tsv", sep='\t')
    df1 = df[ df["Lineage"] == lineage ]
    
    data={}
    labels=[]
    colors={}
     
    the=1.5
    mo=datetime.now().month
    y=datetime.now().year
    mo=5
    y=2021
    
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
            
            
    co=0
    for m in range(mo-5, mo+1):
        df=df1[ ( (df1['Month']==m) & (df1['Year']==y) ) ]
        if(len(df)>0):
            df=df.sort_values('Entropy', ascending=False)
            lab=str(m)+"/"+str(y)
                        
            data[lab]=[ [], [], { 'color': back[co] } ]
            
            pos = []
            end = len(df)
            #if(end>20):
            #    end=20
                
            for i in range(len(df)):
                y_=df.iloc[i,8] # for count
                y_=df.iloc[i,5] # for entropy
                if(y_ > the):
                    p=str(df.iloc[i,4])
                    if(not p in pos):
                        prot=df.iloc[i,3]
                        xlab=prot+"_"+p
                        #if(xlab in cnt.keys()):
                        #    if(cnt[xlab]>2):
                        pos.append(p)
                        
                        data[lab][0].append(xlab)
                        data[lab][1].append(str(y_))
                
                if(len(pos)==3):
                    break
        co+=1                
                    
    bars=[]
    for k in data.keys():
        bars.append( { "x": data[k][0], "y": data[k][1], "name": k, "marker": data[k][2], "type": "bar" } )
     
    resp["error"]=""
    resp["data"]=bars
    
    return json.dumps(resp)

@app.route('/entropy_position_analysis_by_period/<lineage>/<protein>/<position>/<effect>/<month>/<year>', methods=['GET'])
def section_entropy_position_analysis_by_period(lineage, protein, position, effect, month, year):
    groups={ 'Non Polar': ['G','A','V','C','P','L', 'I', 'M','W','F'], 'Polar': ['S','T','Y','N','Q'], 'Positive Charge': ['K','R','H'], 'Negative Charge': ['D','E'] }
    revg={}
    for g in groups:
        for aa in groups[g]:
            revg[aa]=g
            
    uui=str(uuid.uuid4())
    
    data={"error": "", "table_entropy": [], "file_export": "/static/results/"}
    
    df=pd.read_csv(directory+'data/report_position_probabilities_by_lineage_period.tsv', sep='\t')
    
    # Geenerating table
    sr={}
    with open(directory+'data/snap_data.json', 'r') as fp:
        sr = json.load(fp) 
        
    conditions=[]
    lineage=lineage.replace("-",".")
    if(lineage!="All"):
        conditions.append(' ( df["Lineage"]=="'+lineage+'" ) ')
    
    if(protein!="All"):
        cproteins=[]
        for m in protein.split(','):
            cproteins.append(' ( df["Protein"]=="'+m+'" ) ')
        conditions.append("( "+' | '.join(cproteins)+" )")
        
    if(position!="All"):
        conditions.append(' ( df["Position"]=='+position+' ) ')
    
    if(month!="All"):
        conditions.append(' ( df["Month"]=='+month+'  )')
    
    if(year!="All"):
        conditions.append(' ( df["Year"]=='+year+' ) ')
        
    condition="( "+' & '.join(conditions)+" )"
    if(len(conditions)>0):
        df = df[ eval(condition) ]
        
    if(len(df)!=0):
        g=open(directory+'data_exportation/'+uui+".tsv", "w")
        columns=[]
        data['columns']=[]
        for c in df.columns:
            columns.append(c)
            data['columns'].append({ "title" : c } )
        
        columns+=["SNAP effect","SNAP prediction"]
        data['columns']+=[{ "title" : 'SNAP effect' }, { "title": 'SNAP prediction' }]
        
        g.write( ('\t'.join(columns))+"\n")
        
        for i in range(len(df)):
            protein=df.iloc[i,3].lower()
            position = str(df.iloc[i,4])
            aas = df.iloc[i,7].split(" to ")
            
            sn=["-","-"]
            if( protein+"_"+aas[0]+position+aas[1] in sr.keys()):
                sn=sr[protein+"_"+aas[0]+position+aas[1]]
            
            if(sn[0]==effect or effect=="All"):
                dat=list(df.iloc[i, :])+sn
                if(i<1000):
                    data["table_entropy"].append( dat )
                
                c=0
                for d in dat:
                    dat[c]=str(d)
                    c+=1 
                g.write( ('\t'.join(dat))+"\n")
        g.close()
        
        data['file_export'] += uui+'.tsv'
    else:
        data['error']='There are no matches for your filters.'
        
    return json.dumps(data)

# Page structural & functional analysis
@app.route('/structural_functional_analysis', methods=['POST'])
def section_structural_analysis():
    groups={ 'Non Polar': ['G','A','V','C','P','L', 'I', 'M','W','F'], 'Polar': ['S','T','Y','N','Q'], 'Positive Charge': ['K','R','H'], 'Negative Charge': ['D','E'] }
    revg={}
    for g in groups:
        for aa in groups[g]:
            revg[aa]=g
            
    uui=str(uuid.uuid4())
    
    seq = request.form.get('sequence')
    user_email = request.form.get('email')
    
    f=open(directory+"structural_effects/data/ref_sequence.fasta","r")
    for line in f:
        if(line.find(">")==-1):
            ref=line.replace("\n","")
    f.close()
    
    flag=True
    if(seq.find("\n")!=-1):
        st=""
        temp=seq.split("\n")
        for t in temp:
            if(t.find(">")!=-1):
                if(st!=""):
                    flag=False
            
            if(flag and t.find(">")==-1):
                st+=t.replace("*", "")
        seq=st
    simi = Levenshtein.ratio(seq, ref)    
    
    data={"error": "", "error_email": "", "mutations": "", "domains":"", "rmsd": "","table_mutations": [],"table_stability": [], "file_export": "/static/results/"}
    
    if(simi > 0.4):
        os.system("bash "+directory+"structural_effects/run.sh "+seq+" "+uui+" "+directory+"")
        
        sr={}
        snap=pd.read_csv(directory+"structural_effects/spike.csv", sep=",")
        var=list(snap.iloc[:,0])
        effect=list(snap.iloc[:,1])
        accuracy=list(snap.iloc[:,2])
        c=0
        for v in var:
            sr[v]=[effect[c], accuracy[c]]
            c+=1
            
        
        if(os.path.isfile(directory+"structural_effects/"+uui+"/step1_report.txt")):
            f=open(directory+"structural_effects/"+uui+"/step1_report.txt","r")
            for line in f:
                l=line.replace("\n","")
                if(l.find("RMSD")!=-1):
                    data['rmsd']=l.split(":")[1]
            f.close()
            
        if(os.path.isfile(directory+"structural_effects/"+uui+"/mutations.txt")):
            g=open(directory+"structural_effects/"+uui+"/functional_mutations_report.tsv","w")
            g.write("Mutation\tSNAP effect\tSNAP prediction\tAA group reference\tAA group alternative\n")
            
            muts=[]
            pos=[]
            interest=[]
            f=open(directory+"structural_effects/"+uui+"/mutations.txt","r")
            for line in f:
                l=line.replace("\n","")
                
                muts.append(l[1:])
                pos.append(l[2:-1])
                
                ref=revg[l[1]]
                alt="del"
                if(l[1:].find("-")==-1):
                    alt=revg[l[-1]]
                else:
                    sr[l[1:]] = ["-", "-"]
                
                if(sr[l[1:]][0]=='effect'):
                    interest.append(l[1:])
                    
                dat=[ l[1:], sr[l[1:]][0], str(sr[l[1:]][1]), ref, alt ]
                
                data["table_mutations"].append( dat )
                
                g.write(('\t'.join(dat))+"\n")
            f.close()
            
            if(len(pos) > 1):
                data["mutations"] = "("+(' or '.join(pos))+")"
            if(len(pos) == 1):
                data["mutations"] = ' or '.join(pos)
            
            infodom = _get_domain(seq, "Spike") 
            data['domains'] = infodom[0]
            data['domains'] = "<p> <b>Mutation(s) of interest according to functional impact prediction: </b> "+(''.join(interest))+" </p> <br />"
            domsspike=['PS51921-BCOV_S1_CTD','PS51922-BCOV_S1_NTD','PS51923-COV_S2_HR1','PS51924-COV_S2_HR2']
            refcds=['9-303','334-527','896-1001','1143-1225']
            dnotnew=[]
            cnotnew={}
            dnot=[]
            cnot={}
            dok=[]
            cd=0
            for dm in infodom[1]:
                if(dm in domsspike):
                    dok.append(dm+" ("+infodom[2][cd]+")")
                else:
                    dnotnew.append(dm+" ("+infodom[2][cd]+")")
                    if(not dm in cnotnew.keys()):
                        cnotnew[dm]=[dm+" ("+infodom[2][cd]+")", []]
                    cnotnew[dm][1] += _get_muts_in_domains(seq, infodom[2][cd], pos, muts)
                cd+=1
            
            if(len(dok)<len(domsspike)):
                cd = 0
                for dm in domsspike:
                    if(not dm in infodom[1]):
                        dnot.append(dm)
                        if(not dm in cnot.keys()):
                            cnot[dm]=[dm+" ("+refcds[cd]+")", []]
                        cnot[dm][1] += _get_muts_in_domains(seq, refcds[cd], pos, muts)
                    cd+=1
                        
            if(len(dok)>0):
                data['domains']+="<p> <b>Spike Domains in this sequence:</b> "+('; '.join(dok))+"</p> "
                
            if(len(dnot)>0):
                data['domains']+="<p> <b>Spike reference domains not found in this sequence:</b>  </p> <ul>"
                for k in cnot.keys():
                    data['domains']+="<li>"+cnot[k][0]
                    if(len(cnot[k])>0):
                        data['domains']+=" - Mutation(s): "+(', '.join(cnot[k][1]))
                    data['domains']+="</li>"
                data['domains']+="</ul>"
                
            if(len(dnotnew)>0):
                data['domains']+="<p> <b>Extra Domains (not found in Spike reference protein):</b> </p> <ul>"
                for k in cnotnew.keys():
                    data['domains']+="<li>"+cnotnew[k][0]
                    if(len(cnotnew[k])>0):
                        data['domains']+=" - Mutation(s): "+(', '.join(cnotnew[k][1]))
                    data['domains']+="</li>"
                data['domains']+="</ul>"
            
            g.close()
            
            stab={ "BackHbond": "Backbone H. Bonds", "SideHbond": "Sidechain H. Bonds", "Energy_VdW": "Energy Van der Walls", "Electro": "Eletrostatic", "Energy_SolvP": "Energy Solvent Polar", "Energy_SolvH": "Energy Solvent Hydrophobic", "Energy_vdwclash": "Energy Van der Walls clashes", "energy_torsion": "Energy of Torsions", "backbone_vdwclash": "Backbone Van der Walls Clashes", "Entropy_sidec": "Entropy Sidechain", "Entropy_mainc": "Entropy Mainchain", "water bonds": "Water Bonds", "helix dipole": "Helix Dipole", "loop_entropy": "Loop Entropy", "cis_bond": "Cis Bonds", "disulfide": "Disulfide", "kn electrostatic": "Electrostatic Kon", "partial covalent interactions": "Partial Covalent Interactions", "Energy_Ionisation": "Energy of Ionisation", "Entropy Complex": "Entropy of Complex", "Total": "Total Energy"  }
            
            if(os.path.isfile(directory+"structural_effects/"+uui+"/step1_table_stability_comparison.tsv")):
                df=pd.read_csv(directory+"structural_effects/"+uui+"/step1_table_stability_comparison.tsv", sep="\t")
                for i in range(len(df)):
                    
                    if(df.iloc[i, 1]!="backbone_vdwclash" and df.iloc[i,1]!="kn electrostatic"):
                        df.iloc[i, 1] = stab[df.iloc[i, 1][:-1]]
                    else:
                        df.iloc[i, 1] = stab[df.iloc[i, 1]]
                        
                    if(df.iloc[i, 1]=="Total Energy"):
                        diff = df.iloc[i,2]-df.iloc[i,3]
                        if(abs(diff)>50):
                            data['domains']+="<p> <b>According to the energy parameters of stability analysis, the protein has a difference of "+str(diff)+" in relation to the reference Spike protein, which indicates that the structural properties have changed like the solvent accessibility and secondary structure.</b> The complete table of energy stability analysis was sent by e-mail. </p>"
                        else:
                            data['domains']+="<p> <b>According to the energy parameters of stability analysis, the protein has an acceptable difference of total energy in relation to the Spike reference protein, which does not impact in the structural properties.</b> The complete table of energy stability analysis was sent by e-mail.</p>"
                            
                    data["table_stability"].append( list(df.iloc[i, :]) )
            else:
                data["error"] = "The mutations were not found in the reference 3D Spike structure"
        else:
            data["error"] = "There are no mutations in the sequence"
    else:
        data["error"] = "This sequence has less than 40% of similarity with Spike."
            
    shutil.make_archive(directory+"data_exportation/"+uui, 'zip', directory+"structural_effects/"+uui)
    data['file_export'] += uui+'.zip'
    
    os.system("rm -rf "+directory+"structural_effects/"+uui)
    
    try:
        link=host+"/"+data['file_export']
        send_success_email(user_email, link, directory+"data_exportation/"+uui)
    except:
        data["error_email"] = "It was not possible sending the results to the provided e-mail."
        
    return json.dumps(data)

# Page descriptive analysis
@app.route('/get_br_state_plots', methods=['GET'])
def section_get_plots_brstate():
    resp={"error": ""}
    with open(directory+'data/br_state_analysis_plot.json', 'r') as fp:
        dat = json.load(fp)  
    resp["global_geo_info"]=dat
    
    with open(directory+'data/br_state_analysis_data.json', 'r') as fp:
        dat = json.load(fp)  
    resp["global_geo_data"]=dat
    
    return json.dumps(resp)
    
@app.route('/get_plot_forecasting/<location>/<protein>', methods=['GET'])
def section2_get_plot_forecasting(location, protein):
    df=pd.read_csv(directory+'data/forecasting/'+location.replace(" ","-")+"_forecasting.tsv", sep="\t")
    df = df[ df["protein"]==protein ]
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
    plt.savefig(directory+'data_exportation/'+location.replace(" ","-")+"_"+protein+"_forecasting.png")
    
    resp = {"file_export": '/static/results/'+location.replace(" ","-")+"_"+protein+"_forecasting.png"}
    
    return json.dumps(resp)
    
@app.route('/mutations_mean_by_period/<location>/<period>/<protein>/<type_>', methods=['GET'])
def section2_get_mean_mutations_period(location, period, protein, type_):
    # from hive sec0: proteins, locations
    
    # select id_genome_output, nm_nation, nm_lineage, dt_collect, nm_aminoacid_modified, nm_protein from gisaid.vw_genome_protein_domain where nm_nation='Brazil' and nm_protein='Spike'
    
    prs=protein.split(",")
    colors=["#1E90FF","#228B22","#BC8F8F","#B0E0E6"]

    with open(directory+'data/mutation_indel/'+location.replace(" ","-")+"_"+period+'_mutation_indel.json', 'r') as fp:
        dat = json.load(fp)   
    data=dat[type_]    
    labels=dat['labels']
    
    c=0
    #conditions=[]
    for p in prs:
        #conditions.append(" nm_protein='"+p+"' ")

        data[p]["label"]=p
        data[p]["borderColor"]=colors[c]
        data[p]["tension"]=0.1
        data[p]["fill"]=False
        data[p]["backgroundColor"]=colors[c]
        #data[p]["data"]=[]
        #data[p]["deviation"]=[]

        c+=1
    datasets=[]
            
    for p in prs:
        datasets.append(data[p])

    resp={"labels": labels, "datasets": datasets}
    return json.dumps(resp)

@app.route('/lineages_count_by_mutation/<mutations>/<int:hits>', methods=['GET'])
def section3_lineagesCount_by_mutation(mutations, hits): # chart pie
    resp={"error": "Invalid input"}
    try:
        muts=mutations.split(",")
        
        with open(directory+'data/lineageCount_by_mutation.json', 'r') as fp:
            dat = json.load(fp)
            
        counts={}
        for m in muts:
            counts[m] = dat[m]
            
        with open(directory+'data/lineageCount.json', 'r') as fp:
            dat = json.load(fp)
                
        uui=str(uuid.uuid4())
        resp={"error": "", "file_export": "/static/results/"+uui+".tsv"}
        f=open(directory+"data_exportation/"+uui+".tsv","w")
        f.write("Mutation\tLineage\tCount\tProportion\n")
        for m in counts.keys():
            resp[m]={"labels": [], "counts": []}
            sorted_={}
            _list = reversed( sorted( counts[m].items(), key=lambda kv: kv[1] ) )
            i=0
            for l in _list:
                if(i<10 and l[0]!="None"):
                    sorted_[l[0]]=l[1]
                    i+=1
    
            for l in sorted_.keys():
                if(sorted_[l] > hits):
                    resp[m]["labels"].append(l)
                    resp[m]["counts"].append(sorted_[l]/dat[l])
    
            for l in counts[m].keys():
                if(counts[m][l] > hits):
                    f.write("%s\t%s\t%i\t%.2f\n" %( m, l, counts[m][l], counts[m][l]/dat[l] ) )
        f.close()
    except:
        pass

    return json.dumps(resp)

@app.route('/countries_count_by_mutation/<mutations>/<int:hits>', methods=['GET'])
def section5_countriesCount_by_mutation(mutations, hits): # chart pie
    resp={"error": "Invalid input"}
    try:
        muts=mutations.split(",")
    
        with open(directory+'data/countryCount_by_mutation.json', 'r') as fp:
            dat = json.load(fp)
            
        counts={}
        for m in muts:
            counts[m] = dat[m]
    
        with open(directory+'data/locationsCount.json', 'r') as fp:
            dat = json.load(fp)
                
        uui=str(uuid.uuid4())
        resp={"error": "", "file_export": "/static/results/"+uui+".tsv"}
        f=open(directory+"data_exportation/"+uui+".tsv","w")
        f.write("Mutation\tLocation\tCount\tProportion\n")
        for m in counts.keys():
            resp[m]={"labels": [], "counts": []}
            sorted_={}
            _list = reversed( sorted( counts[m].items(), key=lambda kv: kv[1] ) )
            i=0
            for l in _list:
                if(i<10 and l[0]!="None" and l[0] in dat.keys()):
                    sorted_[l[0]]=l[1]
                    i+=1
    
            for l in sorted_.keys():
                if(sorted_[l] > hits):
                    resp[m]["labels"].append(l)
                    resp[m]["counts"].append(sorted_[l]/dat[l])
    
            for l in counts[m].keys():
                if(counts[m][l] > hits and l in dat.keys()):
                    f.write("%s\t%s\t%i\t%.2f\n" %( m, l, counts[m][l], counts[m][l]/dat[l] ) )
        f.close()
    except:
        pass

    return json.dumps(resp)

@app.route('/mutations_count_by_lineage/<lineage>/<int:hits>', methods=['GET'])
def section7_mutationsCount_by_lineage(lineage, hits): # chart pie
    resp={"error": "Invalid input"}
    try:
        muts=lineage.split(",")
        
        with open(directory+'data/mutationCount_by_lineage.json', 'r') as fp:
            dat = json.load(fp)
            
        counts={}
        for m in muts:
            counts[m] = dat[m]
    
        with open(directory+'data/mutationCount.json', 'r') as fp:
            dat = json.load(fp)
                
        uui=str(uuid.uuid4())
        resp={"error": "", "file_export": "/static/results/"+uui+".tsv"}
        f=open(directory+"data_exportation/"+uui+".tsv","w")
        f.write("Lineage\tMutation\tCount\tProportion\n")
        for m in counts.keys():
            resp[m]={"labels": [], "counts": []}
            sorted_={}
            _list = reversed( sorted( counts[m].items(), key=lambda kv: kv[1] ) )
            i=0
            for l in _list:
                if(i<10 and l[0]!="None" and l[0] in dat.keys()):
                    sorted_[l[0]]=l[1]
                    i+=1
    
            for l in sorted_.keys():
                if(sorted_[l] > hits):
                    resp[m]["labels"].append(l)
                    resp[m]["counts"].append(sorted_[l]/dat[l])
    
            for l in counts[m].keys():
                if(counts[m][l] > hits and l in dat.keys()):
                    f.write("%s\t%s\t%i\t%.2f\n" %( m, l, counts[m][l], counts[m][l]/dat[l] ) )
        f.close()
    except:
        pass

    return json.dumps(resp)

@app.route('/unique_mutations_count_by_lineage/<lineage>/<int:hits>', methods=['GET'])
def section8_unique_mutationsCount_by_lineage(lineage, hits): # chart pie
    resp={"error": "Invalid input"}
    try:
        muts=lineage.split(",")
        
        with open(directory+'data/uniqueMutationCount_by_lineage.json', 'r') as fp:
            dat = json.load(fp)
            
        counts={}
        for m in muts:
            counts[m] = dat[m]
    
        with open(directory+'data/mutationCount.json', 'r') as fp:
            dat = json.load(fp)
                
        uui=str(uuid.uuid4())
        resp={"error": "", "file_export": "/static/results/"+uui+".tsv"}
        f=open(directory+"data_exportation/"+uui+".tsv","w")
        f.write("Lineage\tMutation\tCount\tProportion\n")
        for m in counts.keys():
            resp[m]={"labels": [], "counts": []}
            sorted_={}
            _list = reversed( sorted( counts[m].items(), key=lambda kv: kv[1] ) )
            i=0
            for l in _list:
                if(i<10 and l[0]!="None" and l[0] in dat.keys()):
                    sorted_[l[0]]=l[1]
                    i+=1
    
            for l in sorted_.keys():
                if(sorted_[l] > hits):
                    resp[m]["labels"].append(l)
                    resp[m]["counts"].append(sorted_[l]/dat[l])
    
            for l in counts[m].keys():
                if(counts[m][l] > hits and l in dat.keys()):
                    f.write("%s\t%s\t%i\t%.2f\n" %( m, l, counts[m][l], counts[m][l]/dat[l] ) )
        f.close()
    except:
        pass

    return json.dumps(resp)

@app.route('/get_mutations_peptide/<protein>/<peptide>', methods=['GET'])
def section9_get_mutations_peptide(protein, peptide): # chart pie
    resp={"error": "Invalid input", "mutations": ""}
    
    try:
        muts=_get_mutations('test', peptide.upper(), protein)
        
        if(muts==""):
            resp['error']="There is no mutations in the input sequence."
        else:
            with open(directory+'data/mutationCount.json', 'r') as fp:
                dat = json.load(fp)
            
            nodb=[]
            indb=[]
            counts={}
            for m in muts.split(';'):
                if( protein+"_"+m.replace("-", "del") in dat.keys() ):
                    counts[m] = dat[protein+"_"+m.replace("-", "del")]
                    #indb.append(m+" - "+str(counts[m]))
                    indb.append(m)
                else:
                    counts[m] = 0
                    nodb.append(m)
                    
            if(counts=={}):
                resp["error"]="Mutations were not found in database"
            else:
                uui=str(uuid.uuid4())
                resp={"error": "", "file_export": "/static/results/"+uui+".tsv", "labels": [], "counts": []}
                f=open(directory+"data_exportation/"+uui+".tsv","w")
                sorted_={}
                _list = reversed( sorted( counts.items(), key=lambda kv: kv[1] ) )
                i=0
                for l in _list:
                    if(i<10):
                        sorted_[l[0]]=l[1]
                        i+=1
            
                for l in sorted_.keys():
                    resp["labels"].append(l)
                    resp["counts"].append(sorted_[l])
            
                for l in counts.keys():
                    f.write("%s\t%i\n" %( l, counts[l] ) )
                f.close()
                
            resp["mutations_indb"]=indb
            resp["mutations_nodb"]=nodb
    except:
        pass

    return json.dumps(resp)

@app.route('/domain_counts_by_protein', methods=['GET'])
def section6_domain_counts_by_protein(): # chart pie
    fake = Factory.create()
    
    resp={"error": "Invalid input"}
    try:
        with open(directory+'data/info_domains.json', 'r') as fp:
            dat = json.load(fp)
                
        df=pd.read_csv(directory+"data/report_domains.tsv", sep="\t")
        df = df[ df["y"]>50000 ]
        labels=df['label'].unique()
        
        colors=[]
        while (len(colors) != len(labels)):
            color=str(fake.hex_color())
            if(not color in colors):
                colors.append(color)
        
        c=0
        descs={}
        bars = []
        for label, label_df in df.groupby('label'):
            xlab=[]
            for n in label_df.x:
                domain=n.split("-")[1]
                xlab.append(domain)
                descs[domain] = "<p>"+domain+" ("+dat[domain]['prosite']+") - "+dat[domain]['description']+"</p>"
            bars.append( { "x": list(xlab), "y": list(label_df.y), "name": label, "marker": { 'color': colors[c] }, "type": "bar" } )
            c+=1
        
        resp["data"]=bars
        
        resp["error"]=""
        
        resp["info_domain"]=""
        for k in descs.keys():
            resp["info_domain"]+= descs[k]
        
    except:
        pass

    return json.dumps(resp)


if __name__ == "__main__":
    app.run(host=ip)


