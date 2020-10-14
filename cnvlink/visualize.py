# CNVLINK
# CNV2SV project, BCM SV Hackathon October 2020
# M.S.- Mon, October 12 2020

#Prerequisites
#Python >= 3 (Tested on 3.8)
#conda install pandas
#conda install seaborn
#R
#R circlize package

import sys
import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
sns.set()

if len(sys.argv)<2:
    print("usage: visualize.py cnvlink_out_dir")
    exit()
cnvlink_run=sys.argv[1]


def plot_stats(stats,keys,title,plot_name):
    data={"x":[],"y":[]}
    for k,label in keys:
        data["x"].append(label)
        data["y"].append(stats[k]["DUP"])
    #plt.figure(figsize=(6,4))
    fig,ax=plt.subplots(figsize=(6,4))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    sns.barplot(x="x", y="y", data=data, label="Number of CNVs")
    plt.title(title)
    plt.ylabel("Number of CNVs")
    plt.tight_layout()
    plt.savefig("plots/%s_%s.png"%(cnvlink_run,plot_name))
    plt.close()

def plot_dist(data,title,plot_name,y_log=False):
    sns.distplot(data)
    plt.title(title)
    plt.tight_layout()
    plt.savefig("plots/%s_%s.png"%(cnvlink_run,plot_name))
    plt.close()

def plot_box(data,title,plot_name,set_ylim=False,ylabel=None):
    data_df={"source":[],"y":[]}
    for y,source in data:
        for d in y:
            data_df["source"].append(source)
            data_df["y"].append(d)
    df=pd.DataFrame(data_df)
    sns.boxplot(data=df,y="y",x="source",orient="v",showfliers = False)
    sns.swarmplot(data=df,y="y",x="source",orient="v",color="0.25")
    #sns.swarmplot(data=df,y="dist",x="source",orient="v")
    if set_ylim!=False:
        plt.ylim(set_ylim)
    plt.title(title)
    plt.xlabel(None)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig("plots/%s_%s.png"%(cnvlink_run,plot_name))
    plt.close()

import random
def plot_scatter(x,y,title,plot_name,xlabel,ylabel):
    jitter_amnt=0.05
    sns.scatterplot([k+random.uniform(-jitter_amnt,jitter_amnt) for k in x],[k+random.uniform(-jitter_amnt,jitter_amnt) for k in y]) #Jitter
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig("plots/%s_%s.png"%(cnvlink_run,plot_name))
    plt.close()

stats=json.load(open(cnvlink_run+"/cnvlink_stats.json","r"))
stats["adj_total"]={"DUP":stats["matched_adj"]["DUP"]+stats["perfect_adj"]["DUP"]}
stats["far_total"]={"DUP":stats["matched_far"]["DUP"]+stats["perfect_far"]["DUP"]}
stats["link_total"]={"DUP":stats["adj_total"]["DUP"]+stats["far_total"]["DUP"]}
stats["unlinked"]={"DUP":stats["cnv_call"]["DUP"]-stats["adj_total"]["DUP"]-stats["far_total"]["DUP"]}
stats["far_same_chr"]={"DUP":0}
stats["far_other_chr"]={"DUP":0}

far_same_chr_startpos_dists=[]
adj_startpos_dists=[]

far_length_diffs=[]
adj_length_diffs=[]

far_aln_score_frac=[]
adj_aln_score_frac=[]

good_matches_adj_list=[]
good_matches_far_list=[]

for line_ in open(cnvlink_run+"/cnvlink_out.tsv","r").readlines()[3:]:
    line=line_.strip()
    cnv_id,status,fail_reason,type,good_matches_adj,good_matches_far,alignment_matches,cnv_chr,cnv_start,cnv_length,asm_chr,asm_start,asm_length,all_good_matches_adj,all_good_matches_far,aln_cigar,aln_length,aln_NM,aln_mapq,cnv_seq,asm_seq=line.split("\t")

    good_matches_adj_list.append(int(good_matches_adj))
    good_matches_far_list.append(int(good_matches_far))
    if status=="match_far" or status=="perfect_far":
        if cnv_chr==asm_chr:
            stats["far_same_chr"]["DUP"]+=1
            far_same_chr_startpos_dists.append(abs(int(cnv_start)-int(asm_start)))
        else:
            stats["far_other_chr"]["DUP"]+=1
        far_aln_score_frac.append(float(alignment_matches)/float(cnv_length))
        far_length_diffs.append(abs(int(cnv_length)-int(asm_length)))
    if status=="match_adj" or status=="perfect_adj":
        adj_startpos_dists.append(abs(int(cnv_start)-int(asm_start)))
        adj_aln_score_frac.append(float(alignment_matches)/float(cnv_length))
        adj_length_diffs.append(abs(int(cnv_length)-int(asm_length)))

plot_box([(adj_length_diffs,"Adjacent SV"),(far_length_diffs,"Distant SV")],"Length difference between called CNV and linked SV","lengthdiff_box",set_ylim=(0,5500),ylabel="abs(Length difference)")
plot_box([(adj_aln_score_frac,"Adjacent SV"),(far_aln_score_frac,"Distant SV")],"Alignment concordance between called CNV and linked SV","aln_box",ylabel="Matching bp / Length of CNV")
plot_scatter(good_matches_adj_list,good_matches_far_list,"Number of adjacent vs distant SV matches for each CNV","matches_adj_vs_far",xlabel="Number of adjacent linked SVs",ylabel="Number of distant linked SVs")

plot_stats(stats,[("unlinked","Not linked"),("link_total","Linked to SV")],title="Linked SVs are identified for a fraction of CNVs called by CNVnator",plot_name="linked") #,("unlinked","Not linked")])
plot_stats(stats,[("adj_total","Adjacent"),("far_same_chr","Distant"),("far_other_chr","Other Chromosome")],title="Location of CNV-call linked SVs",plot_name="location") #,("unlinked","Not linked")])
