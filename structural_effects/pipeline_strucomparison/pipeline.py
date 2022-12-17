#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 16:21:30 2021

@author: yasmmin
"""

import pandas as pd
import os
import json
import re
import statistics as st

import math

import matplotlib.pyplot as plt
import numpy as np
from Bio import PDB
from matplotlib import colors

class StructuralEffect:
    def _get_domain(self, sequence, protein, prefix):
        with open(prefix+"data/train_domains.json", 'r') as fp:
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
            
        return "; ".join(text)

    def plot_ramachandran(self, file):
        __file__=file
    
        """
        The preferences were calculated from the following artice:
        Lovell et al. Structure validation by Calpha geometry: phi,psi and Cbeta deviation. 2003
        DOI: 10.1002/prot.10286
        """
    
        # General variable for the background preferences
        rama_preferences = {
            "General": {
                "file": "rama-500/rama500-general.data",
                "cmap": colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
                "bounds": [0, 0.0005, 0.02, 1],
            },
            "GLY": {
                "file": "rama-500/rama500-gly-sym.data",
                "cmap": colors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F']),
                "bounds": [0, 0.002, 0.02, 1],
            },
            "PRO": {
                "file": "rama-500/rama500-pro.data",
                "cmap": colors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
                "bounds": [0, 0.002, 0.02, 1],
            },
            "PRE-PRO": {
                "file": "rama-500/rama500-prepro.data",
                "cmap": colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
                "bounds": [0, 0.002, 0.02, 1],
            }
        }
    
        # Read in the expected torsion angles
        __location__ = './' #You must set the ptah of the .data files here
        rama_pref_values = {}
        for key, val in rama_preferences.items():
            rama_pref_values[key] = np.full((360, 360), 0, dtype=np.float64)
            with open(os.path.join(__location__, val["file"])) as fn:
                for line in fn:
                    if not line.startswith("#"):
                        # Preference file has values for every second position only
                        rama_pref_values[key][int(float(line.split()[1])) + 180][int(float(line.split()[0])) + 180] = float(
                            line.split()[2])
                        rama_pref_values[key][int(float(line.split()[1])) + 179][int(float(line.split()[0])) + 179] = float(
                            line.split()[2])
                        rama_pref_values[key][int(float(line.split()[1])) + 179][int(float(line.split()[0])) + 180] = float(
                            line.split()[2])
                        rama_pref_values[key][int(float(line.split()[1])) + 180][int(float(line.split()[0])) + 179] = float(
                            line.split()[2])
    
        normals = {}
        outliers = {}
        for key, val in rama_preferences.items():
            normals[key] = {"x": [], "y": [],'Res':[]}
            outliers[key] = {"x": [], "y": []}
    
        # Calculate the torsion angle of the inputs
        structure = PDB.PDBParser().get_structure('input_structure', __file__)
        for model in structure:
            for chain in model:
                polypeptides = PDB.PPBuilder().build_peptides(chain)
                for poly_index, poly in enumerate(polypeptides):
                    phi_psi = poly.get_phi_psi_list()
                    for res_index, residue in enumerate(poly):
                        res_name = "{}".format(residue.resname)
                        res_num = residue.id[1]
                        phi, psi = phi_psi[res_index]
                        if phi and psi:
                            aa_type = ""
                            if str(poly[res_index + 1].resname) == "PRO":
                                aa_type = "PRE-PRO"
                            elif res_name == "PRO":
                                aa_type = "PRO"
                            elif res_name == "GLY":
                                aa_type = "GLY"
                            else:
                                aa_type = "General"
                            if rama_pref_values[aa_type][int(math.degrees(psi)) + 180][int(math.degrees(phi)) + 180] < \
                                    rama_preferences[aa_type]["bounds"][1]:
                                print("{} {} {} {}{} is an outlier".format(inp, model, chain, res_name, res_num))
                                outliers[aa_type]["x"].append(math.degrees(phi))
                                outliers[aa_type]["y"].append(math.degrees(psi))
                            else:
                                normals[aa_type]["x"].append(math.degrees(phi))
                                normals[aa_type]["y"].append(math.degrees(psi))
                                normals[aa_type]['Res'].append(res_name+'_'+str(res_num))
    
        # Generate the plots
        plt.figure(figsize=(10,10))
        for idx, (key, val) in enumerate(sorted(rama_preferences.items(), key=lambda x: x[0].lower())):
            plt.subplot(2, 2, idx + 1)
            plt.title(key,fontsize=20)
            plt.imshow(rama_pref_values[key], cmap=rama_preferences[key]["cmap"],
                       norm=colors.BoundaryNorm(rama_preferences[key]["bounds"], rama_preferences[key]["cmap"].N),
                       extent=(-180, 180, 180, -180),alpha=0.7)
    
            plt.scatter(normals[key]["x"], normals[key]["y"],s=[15],marker='.')
    
            #for key in normals:
                #for i, name in enumerate (normals[key]['Res']):
                    #plt.annotate(name, (normals[key]["x"][i], normals[key]["y"][i]))
    
            plt.scatter(outliers[key]["x"], outliers[key]["y"],color="red",s=[15],marker='.')
            plt.xlim([-180, 180])
            plt.ylim([-180, 180])
            plt.plot([-180, 180], [0, 0],color="k",alpha=0.7)
            plt.plot([0, 0], [-180, 180],color="k",alpha=0.7)
            plt.xlabel(r'$\phi$',fontsize=12)
            plt.ylabel(r'$\psi$',fontsize=12)
            plt.grid(linestyle='dotted')
    
        plt.tight_layout()
        name=__file__.replace(".pdb","")+"_plot_ramachandran"
        plt.savefig(name+".png", dpi=300)
        
    def _get_mutations(self, id_, seq, uui, prefix):
        f=open(uui+".fasta","w")
        f.write(">"+str(id_)+"\n"+seq+"\n")
        f.close()
        
        mutations=""
        if(os.path.isfile(prefix+"structural_effects/pipeline_strucomparison/blastdb/spike_ref.phr")):
            os.system("blastp -db "+prefix+"structural_effects/pipeline_strucomparison/blastdb/spike_ref -query "+prefix+"structural_effects/pipeline_strucomparison/"+uui+".fasta -num_alignments 1 -out "+prefix+"structural_effects/pipeline_strucomparison/"+uui+".out -outfmt 3")
            
            refseq=""
            stref=0
            
            qseq=""
            stqseq=1
            
            found=False
            init=0
            co=0
            
            changes=[]
            f=open(prefix+"structural_effects/pipeline_strucomparison/"+uui+".out","r")
            for line in f:
                if(line.startswith("Length")):
                    aux=line.replace("\n","").split("=")[1]
                    init=len(aux)+11
                    
                if(line.startswith("0 ") and not found):
                    l=line.replace("\n","").split(" ")
                    elements=[]
                    for el in l:
                        if(el!=""):
                            elements.append(el)
                    stref=int(elements[1])
                    refseq=line.replace("\n","")[init:init+60]
                    stqseq=0
                    for s in qseq:
                        if(s==' '):
                            break
                        else:
                            if(s!='X' and s!='N' and s!='n'):
                                if(refseq[stqseq] != "." and refseq[stqseq] != " " and (s!='-' or refseq[stqseq] != "-") ):
                                    changes.append( "A"+(refseq[stqseq])+str(stref)+s )
                            stref+=1
                            stqseq+=1
                            
                    found=True
                    
                if(line.startswith("Query_1")):
                    found=False
                    qseq=line.replace("\n","")[init:init+60]
                    
                co+=1
            f.close()
            
            mutations=';'.join(changes)
            
        os.system("rm "+prefix+"/structural_effects/pipeline_strucomparison/"+uui+".*")
        
        return mutations
        
    def step1_generate_mutated_structure(self, folder, ref_pdb, list_mutations, seq, prefix): # input ref_sequence, ref_pdb, list_mutations
        print("Step 1 - Generation of mutated structure")
        
        uui=folder.split("/")[-2]
        
        if(list_mutations=='None'):
            muts=self._get_mutations('test', seq, uui, prefix)
            list_mutations=""
            
            if( len(muts.split(";")) > 0 and muts!=""):
                f=open(folder+"mutations.txt","w")
                f.write("\n".join(muts.split(";"))+"\n")
                f.close()
                list_mutations="mutations.txt"
        
        found = {}
        if(list_mutations!=""):
            data=[]
            f=open(folder+ref_pdb, "r")
            for line in f:
                features=[]
                if(line.startswith("ATOM")):
                    l=line.replace("\n","").split(" ")
                    for el in l:
                        if(el!=""):
                            features.append(el)
                    data.append(features)
            f.close()
            
            mutations=[]
            complete=[]
            f=open(folder+list_mutations, "r")
            for line in f:
                l=line.replace("\n","")
                
                mutations.append(l[:-1])
                complete.append(l[1]+l[0]+l[2:])
            f.close()
            
            d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
            
            df = pd.DataFrame(data=data)
            
            for i in range(len(df)):
                try:
                    chain=df.iloc[i,4]
                    ref=d[df.iloc[i,3]]
                    pos=df.iloc[i,5]
                    value=chain+ref+pos
                    if( (value in mutations) ):
                        found[value]=0
                except:
                    pass
            
        found=list(found.keys())  
        
        g=open(folder+"step1_report.txt","w")
        if(len(found)==0):
            g.write("Error: The reference structure do not have the positions of the required mutations.\n")
            g.write("Step 1 will not continue.\n")
        else:
            if( os.path.isdir(folder+"step1_out")):
                os.system("rm -rf "+folder+"step1_out")
            os.system("mkdir "+folder+"step1_out")
            
            diff=len(mutations)-len(found)
            ok=[]
            no=[]
            c=0
            for m in mutations:
                if(m in found ):
                    if(not complete[c] in ok):
                        ok.append(complete[c])
                else:
                    if(not complete[c] in no):
                        no.append(complete[c])
                c+=1
            g.write("Mutations in reference structure: %s\n" %(','.join(ok)) )
            if(len(no)>0):
                g.write("Mutations not found in reference structure: %s\n" %(','.join(no)) )
                g.write("Warning: Step 1 will continue with the valid required mutations\n")
            
               
            # Writing mutation file for foldx
            f=open(folder+"step1_out/individual_list.txt","w")
            f.write("%s;" %(";".join(ok)) )
            f.close()
            
            alternative = ref_pdb.replace(".pdb","")+"_1.pdb"
            print("    Generating mutated structure")
            os.system(prefix+"structural_effects/pipeline_strucomparison/foldx/foldx --command BuildModel --pdb "+ref_pdb+" --pdb-dir "+folder+" --output-dir "+folder+" --output-file step1_modelling_out --mutant-file "+folder+"step1_out/individual_list.txt")
            os.system("mv "+folder+"*.fxout "+folder+"step1_out")
            
            g.write("    Mutated 3D structure saved in "+folder+alternative)
            end="END                                                                             "
            new_pdb=[]
            f=open(folder+alternative,"r")
            for line in f:
                if(line.startswith("ATOM")):
                    new_pdb.append(line.strip())
            f.close()
            
            f=open(folder+alternative,"w")
            f.write("\n".join(new_pdb))
            f.write("\n"+end)
            f.close()
            
            print("Generating Ramachandran plot")
            self.plot_ramachandran(folder+alternative)
            
            g.write("Results for structure "+alternative+":\n")
            
            print("    Calculating RMSD")
            os.system(prefix+"structural_effects/pipeline_strucomparison/foldx/foldx --command Rmsd --pdb-dir "+folder+" --pdb1 "+ref_pdb+" --pdb2 "+alternative+" > "+folder+"step1_out/step1_rmsd_out ")
            rmsd=""
            f=open(folder+"step1_out/step1_rmsd_out","r")
            for line in f:
                if(line.find("rmsd:")!=-1):
                    l=line.replace("rmsd:", "").split(" ")
                    for el in l:
                        if(el!=""):
                            rmsd=el
            f.close()
            g.write("    RMSD: "+rmsd+"\n")
            
            domains = self._get_domain(seq, "Spike", prefix)
            g.write("    Found domains: "+domains+"\n")
            
            print("    Computing stability")
            stability={}
            os.system(prefix+"structural_effects/pipeline_strucomparison/foldx/foldx --command Stability --pdb-dir "+folder+" --pdb "+ref_pdb+" > "+folder+"step1_out/step1_stability_reference_out ")
            f=open(folder+"step1_out/step1_stability_reference_out","r")
            for line in f:
                if(line.find("=")!=-1):
                    l=line.replace("\n", "").split("=")
                    stability[re.sub(' +', ' ', l[0])] = [ float(l[1].replace(" ","")) ]
            f.close()
            
            os.system(prefix+"structural_effects/pipeline_strucomparison/foldx/foldx --command Stability --pdb-dir "+folder+" --pdb "+alternative+" > "+folder+"step1_out/step1_stability_mutated_out ")
            f=open(folder+"step1_out/step1_stability_mutated_out","r")
            for line in f:
                if(line.find("=")!=-1):
                    l=line.replace("\n", "").split("=")
                    stability[re.sub(' +', ' ', l[0])].append( float(l[1].replace(" ","")) )
            f.close()
            
            g.write("    Stability: see table in step1_table_stability_comparison.tsv\n")
            gf=open(folder+"step1_table_stability_comparison.tsv", "w")
            gf.write("Decision\tVariable\tReference\tMutated\n")
            for k in stability.keys():
                decision="same"
                if( len(stability[k]) > 1):
                    if(stability[k][0] > stability[k][1]):
                        decision = "decreased"
                    if(stability[k][0] < stability[k][1]):
                        decision = "increased"
                    
                    gf.write("%s\t%s\t%.3f\t%.3f\n" %(decision, k, stability[k][0], stability[k][1]) )
            gf.close()   
            
            os.system("rm ./*.fxout")
        g.close()
        
    def step2_binding_pockets(self, folder, ref_pdb):
        print("Step 2 - Calculating binding pockets")
        
        alternative_pdb = ref_pdb.replace(".pdb","")+"_1.pdb"
        if(os.path.isfile(folder+alternative_pdb)):
            g=open(folder+"step2_report.txt","w")
            g.write("Binding sites and druggable pockets analysis\n")
            
            if( os.path.isdir(folder+"step2_out")):
                os.system("rm -rf "+folder+"step2_out")
            os.system("mkdir "+folder+"step2_out")
            
            dir_ref=ref_pdb.split(".")[0]+"_out"
            pockets=ref_pdb.split(".")[0]+"_info.txt"
            os.system("./fpocket/bin/fpocket -f "+folder+ref_pdb)
            os.system("mv "+folder+dir_ref+" "+folder+"step2_out/")
            
            i=0
            id_=""
            ref={}
            druggable={}
            f=open(folder+"step2_out/"+dir_ref+"/"+pockets, "r")
            for line in f:
                line=line.replace("\n","")
                if(line.find("    ")!=-1):
                    l=line.replace(":","").split("    ")
                    ref[id_][re.sub(' +', ' ', l[1])]=l[2]
                    if(l[1].find("Druggability Score")!=-1):
                        if( float(l[2]) > 0.7):
                            x=[]
                            y=[]
                            z=[]
                            gf=open(folder+"step2_out/"+dir_ref+"/pockets/pocket"+id_.replace("P","")+"_atm.pdb","r")
                            for line2 in gf:
                                if(line2.find("ATOM")!=-1):
                                    fea=[]
                                    l2=line2.replace("\n","").split(" ")
                                    for el in l2:
                                        if(el!=""):
                                            fea.append(el)
                                    x.append(float(fea[6]))
                                    y.append(float(fea[7]))
                                    z.append(float(fea[8]))
                            gf.close()
                            
                            x=str(st.mean(x))
                            y=str(st.mean(y))
                            z=str(st.mean(z))
                            
                            druggable[id_]=[x,y,z]
                else:
                    if(line!=""):
                        id_="P"+str(i+1)
                        ref[id_]={}
                        i+=1
            f.close()
            
            g.write("    Results for reference structure:\n")
            g.write("        Druggable pockets: see their center coordinates in step2_table_druggable_reference.tsv\n")
            gf=open(folder+"step2_table_druggable_reference.tsv", "w")
            gf.write("Id    center x    center y    center z\n")
            for k in druggable.keys():
                gf.write("%s    %s\n" %(k, '    '.join(druggable[k]) ) )
            gf.close()
            
            g.write("        All pockets: see their center coordinates in step2_table_pockets_reference.tsv\n")
            gf=open(folder+"step2_table_pockets_reference.tsv", "w")
            gf.write("Id    %s\n" %('    '.join(list(ref["P1"].keys())) ) )
            for k in ref.keys():
                gf.write("%s    %s\n" %(k, '    '.join(ref[k].values()) ) )
            gf.close()
            
            dir_ref=alternative_pdb.split(".")[0]+"_out"
            pockets=alternative_pdb.split(".")[0]+"_info.txt"
            os.system("./fpocket/bin/fpocket -f "+folder+alternative_pdb)
            os.system("mv "+folder+dir_ref+" "+folder+"step2_out/")
            
            i=0
            id_=""
            ref={}
            druggable={}
            f=open(folder+"step2_out/"+dir_ref+"/"+pockets, "r")
            for line in f:
                line=line.replace("\n","")
                if(line.find("    ")!=-1):
                    l=line.replace(":","").split("    ")
                    ref[id_][re.sub(' +', ' ', l[1])]=l[2]
                    if(l[1].find("Druggability Score")!=-1):
                        if( float(l[2]) > 0.7):
                            x=[]
                            y=[]
                            z=[]
                            gf=open(folder+"step2_out/"+dir_ref+"/pockets/pocket"+id_.replace("P","")+"_atm.pdb","r")
                            for line2 in gf:
                                if(line2.find("ATOM")!=-1):
                                    fea=[]
                                    l2=line2.replace("\n","").split(" ")
                                    for el in l2:
                                        if(el!=""):
                                            fea.append(el)
                                    x.append(float(fea[6]))
                                    y.append(float(fea[7]))
                                    z.append(float(fea[8]))
                            gf.close()
                            
                            x=str(st.mean(x))
                            y=str(st.mean(y))
                            z=str(st.mean(z))
                            
                            druggable[id_]=[x,y,z]
                else:
                    if(line!=""):
                        id_="P"+str(i+1)
                        ref[id_]={}
                        i+=1
            f.close()
            
            g.write("    Results for mutated structure:\n")
            g.write("        Druggable pockets: see their center coordinates in step2_table_druggable_mutated.tsv\n")
            gf=open(folder+"step2_table_druggable_mutated.tsv", "w")
            gf.write("Id    center x    center y    center z\n")
            for k in druggable.keys():
                gf.write("%s    %s\n" %(k, '    '.join(druggable[k]) ) )
            gf.close()
            
            g.write("        All pockets: see their center coordinates in step2_table_pockets_mutated.tsv\n")
            gf=open(folder+"step2_table_pockets_mutated.tsv", "w")
            gf.write("Id    %s\n" %('    '.join(list(ref["P1"].keys())) ) )
            for k in ref.keys():
                gf.write("%s    %s\n" %(k, '    '.join(ref[k].values()) ) )
            gf.close()
            
            g.close()
        else:
            print("Error: the mutated structure was not found in the working directory")
            
    def step3_secondary_structure_features(self, folder, ref_sequence, list_mutations):
        print("Step 3 - Calculating divergence on secondary structure features")
        
        try:
            g=open(folder+"step3_report.txt","w")
            g.write("    Generating mutated amino acid sequence\n")
            g.close()
            
            ref=""
            f=open(folder+ref_sequence, "r")
            for line in f:
                l=line.replace("\n","")
                if(not l.startswith(">")):
                    ref+=l
            f.close()
            
            mutated_sequence=list(ref)
            f=open(folder+list_mutations, "r")
            for line in f:
                l=line.replace("\n","")
                pos=int(l[2:-1])
                change=l[-1]
                
                mutated_sequence[pos-1]=change
            f.close()
            mutated_sequence = ''.join(mutated_sequence)
            
            if( os.path.isdir(folder+"step3_out")):
                os.system("rm -rf "+folder+"step3_out")
            os.system("mkdir "+folder+"step3_out")
            
            gf=open(folder+"step3_out/mutated_sequence.fasta","w")
            gf.write(">mutated\n%s\n" %(mutated_sequence) )
            gf.close()
            mutated_sequence=folder+"step3_out/mutated_sequence.fasta"
            
            with open(folder+"step3_report.txt","a") as g:
                g.write("    Predicting 2D structure features for reference sequence\n")
            name_ref=ref_sequence.split(".")[0]
            folder_ref=folder+"step3_out/reference_properties"
            os.system("./Predict_Property/Predict_Property.sh -i "+folder+ref_sequence+" -o "+folder+"step3_out/reference_properties")
            
            with open(folder+"step3_report.txt","a") as g:
                g.write("    Predicting 2D structure features for mutated sequence\n")
            name_mut='mutated_sequence'
            folder_mut=folder+"step3_out/mutated_properties"
            os.system("./Predict_Property/Predict_Property.sh -i "+folder+mutated_sequence+" -o "+folder+"step3_out/mutated_properties")
            
            names=["Types of 2D structures", "Disorder","Solvent accessibility"]
            features=["ss3","diso","acc"]
            c=0
            for fea in features:
                with open(folder+"step3_report.txt","a") as g:
                    g.write("    Comparison results for feature "+names[c]+": see table step3_table_"+fea+"_comparison.tsv\n")
                self.auxStep3_parse_2dfeatures( folder, folder_ref, name_ref, folder_mut, name_mut, fea)
                c+=1
                
        except:
            print("Error: the specified mutation does not match the provided reference sequence")
        
    def auxStep3_parse_2dfeatures(self, folder, folder_ref, name_ref, folder_mut, name_mut, feature):
        meanings={ "ss3": { "H": "Helix", "E": "Strand", "C": "Coil" }, "diso": { "*": "disordered", ".": "ordered" }, "acc": { "B": "Bury", "M": "Medium", "E": "Exposed" } }
        
        count={"ref": {}, "mut": {} }
        features={"ref": [], "mut": [] }
        percent={"ref": [], "mut": [] }
        
        strs=["Reference", "Mutated"]
        co=0
        for k in count.keys():
            c=0
            f=open(eval("folder_"+k)+"/"+eval("name_"+k)+"."+feature+"_simp")
            for line in f:
                if(c==2):
                    l=line.replace("\n","")
                    for el in l:
                        if(not el in count[k].keys()):
                            count[k][el]=0
                        count[k][el]+=1
                        features[k].append(el)
                c+=1
            f.close()
            
            total=sum(count[k].values())
            for el in count[k]:
                percent[k].append(el+" ("+meanings[feature][el]+") - "+str(count[k][el]/total))
                
            with open(folder+"step3_report.txt","a") as g:
                g.write("        %s: %s\n" %(strs[co], ', '.join(percent[k]) ) )
                    
            co+=1
            
        c=0
        f=open(folder+"step3_table_"+feature+"_comparison.tsv", "w")
        f.write("decision    reference    mutated\n")
        for k in features["ref"]:
            decision="same"
            if(features["mut"][c] != k):
                decision="changed"
            
            f.write("%s    %s    %s\n" %(decision, k+"-"+meanings[feature][k], features["mut"][c]+"-"+meanings[feature][features["mut"][c]]) )
            c+=1
        f.close()
                
class Running_config:
    def __init__(self):
        self.obj=StructuralEffect()

    def run_step1_mode1(self, folder, ref_pdb, list_mutations, seq, prefix):
        self.obj.step1_generate_mutated_structure(folder, ref_pdb, list_mutations, seq, prefix)

    def run_step2_mode1(self, folder, ref_pdb):
        self.obj.step2_binding_pockets(folder, ref_pdb)

    def run_step3_mode1(self, folder, ref_sequence, list_mutations):
        self.obj.step3_secondary_structure_features(folder, ref_sequence, list_mutations)

    def run(self, args):
        if(args.folder!="" and os.path.isdir(args.folder)):
            #print(args.file_evaluation, args.folder)
            run=0
            if(args.running_type=="" ):
                run=0
            else:
                if(args.running_type in [0,1,2,3]):
                    run=args.running_type
                else:
                    print("Error: invalid choice")
            
            if(args.reference_pdb=="" or args.reference_sequence=="" or args.list_mutations==""):
                print("Error: You have to specify the required files: reference pdb, reference amino acid sequence and list of mutations")
            else:
                folder=args.folder
                ref_pdb=args.reference_pdb
                ref_sequence=args.reference_sequence
                list_mutations=args.list_mutations
                
                if(run==0 or run==1):
                    print("Running step 1")
                    self.run_step1_mode1(folder, ref_pdb, list_mutations, ref_sequence, args.prefix)
                    
                    if(run==0):
                        print("Running step 2")
                        self.run_step2_mode1(folder, ref_pdb)
                        print("Running step 3")
                        self.run_step3_mode1(folder, ref_sequence, list_mutations)
                else:
                    if(run==2):
                        self.run_step2_mode1(folder, ref_pdb)
                    else:
                        self.run_step3_mode1(folder, ref_sequence, list_mutations)

        else:
            print("Error: You have to specify a valid folder to store files")

# call python3 pipeline.py -fo ../data/ -rp nsp5_original.pdb -rf ref_sequence.fasta -lm mutations.txt -rt 1

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description='Structural Features Comparison', formatter_class=RawTextHelpFormatter)
parser.add_argument("-fo", "--folder", action="store", help="(For both modes) Folder to store the files (use the folder where the required files can be found, ex.: /home/user/experiment/ )\n")

parser.add_argument("-rt", "--running_type", action="store", help="0 (default) - Run all steps\n\
1 - Run step 1 (Generate mutated 3D structure)\n\
2 - Run step 2 (Calculate binding sites and possible druggable pockets)\n\
3 - Run step 3 (Analyze 2D structure properties)", type=int)
parser.add_argument("-rp", "--reference_pdb", action="store", help="PDB file of the reference protein")
parser.add_argument("-rf", "--reference_sequence", action="store", help="Fasta file of the reference protein\n")
parser.add_argument("-lm", "--list_mutations", action="store", help="Text file with each mutation by line, ex.: AD614G, where A is the chain, D is the reference amino acid, 614 is the position and G is the changed amino acid")

parser.add_argument("-pr", "--prefix", action="store", help="Prefix\n")

args = parser.parse_args()
r=Running_config()
r.run(args)  
