# CNVLINK
# CNV2SV project, BCM SV Hackathon October 2020
# M.S.- Mon, October 12 2020

#Prerequisites
#Python >= 3 (Tested on 3.8)
#conda install intervaltree
#conda install biopython
#conda install samtools

#TODO: Align unmatched CNV sequences against dup/del sequences from assembly-assembly vcf file (distant matching), this will allow checking for non-adjacent CNVs
#TODO: Confirm results

import os
import sys
import intervaltree
import collections
import subprocess
import json
from Bio import pairwise2
from Bio.Seq import Seq

version="1.0b"

cnv_vcf_file=sys.argv[1]       #CNV calls output vcf file (calls made on GRCH38 with chm13 reads)
asm_asm_vcf_file=sys.argv[2]   #dipcall output vcf file (CHM13 vs GRCH38 assembly-assembly alignment)
ref_file=sys.argv[3]           #Reference file (GRCH38)
output_dir=sys.argv[4].strip() #Output directory to store the resulting files

svtype_filter=["DUP"]

#Parameters for confirming CNVs where duplications / deletions are adjacent
padding_bp=1000           #Padding around CNV calls added to starts/end position when recorded in intervaltree, larger value allows for larger distance between called CNV and asm-asm event
min_cnv_length=50         #Minimum bp for a variant called in asm_asm_vcf_file to be considered as a potential duplication/deletion event
max_cnv_length=1000*50    #Maximum length of CNVs to be included in the analysis (performance)
max_length_diff_frac=0.2  #Maximum length difference when matching called CNV to asm-asm events (in terms of fraction of length of longer sequence)
max_start_pos_diff=1000   #Maximum difference in starting position for adjacent duplication events
min_score_frac=0.8        #Minimum global alignment identity between CNV and asm-asm event sequence

stats_adj={
       "putative":      {"DUP":0,"DEL":0},  #Total number of asm-asm putative dup/del events
       "cnv_call":      {"DUP":0,"DEL":0},  #Total number of CNV calls
       "multi_overlap": {"DUP":0,"DEL":0},  #CNVs that overlap multiple asm-asm events (currently, these CNVs are skipped!)
       "matched":       {"DUP":0,"DEL":0},  #CNVs with a corresponding asm-asm event that fulfills all criteria
       "perfect":       {"DUP":0,"DEL":0},  #CNVs with a corresponding asm-asm event that fulfills all criteria and exactly match in terms of start position and length
       "bad_nonsense":  {"DUP":0,"DEL":0},  #e.g. matching DUP with DEL
       "bad_startdiff": {"DUP":0,"DEL":0},  #Start position difference too large
       "bad_alnscore":  {"DUP":0,"DEL":0},  #Alignment score between CNV sequence and asm-asm event too low
       "bad_lendiff":   {"DUP":0,"DEL":0},  #Event length difference too large
       "bad_no_overlap":{"DUP":0,"DEL":0},  #CNVs with no overlapping asm-asm event
       "bad_other":     {"DUP":0,"DEL":0}   #All other match failure reasons
       }


def readvcf(filename):
    for line_ in open(filename,"r"):
        line=line_.strip()
        if line=="" or line[0]=="#":
            continue
        chr,pos,id,ref,alt,qual,filter,info_str,format,extra=line.split("\t")[:10]
        var={"chr":chr,"pos":int(pos),"id":id,"ref":ref,"alt":alt,"qual":qual,"filter":filter,"info_str":info_str,"format":format,"extra":extra}
        yield var

def vcfinfo2dict(info_str):
    result={}
    for part in info_str.split(";"):
        if "=" in part:
            key,value=part.split("=")
            result[key]=value
    return result

def extractseq(filename,chr,start,end,header=None):
    cmd="samtools faidx %s %s:%d-%d"%(filename,chr,start,end)
    proc=subprocess.run(cmd, shell=True, capture_output=True, check=True)
    if proc.stdout==None:
        return None
    out=proc.stdout.decode("utf8")
    #n_count=sum([1 if c=="N" or c=="n" else 0 for c in out])
    #if float(n_count)/len(out) > max_n_frac:
    #    return None

    seq="".join(out.split("\n")[1:])+"\n"
    if len(seq.strip())=="":
        return None
    if header!=None:
        return ">%s"%(id)+"\n"+seq
    else:
        return seq.strip()


def evalmatch(cnv,asm):
    #Matched e.g. a CNV duplication call with a deletion event
    if asm["putative_svtype"] != cnv["svtype"]:
        #cnv["reason"]="nonsense"
        #stats_adj["nonsense"][cnv["svtype"]]+=1
        return False,"bad_nonsense"

    #Length difference of reported event too long
    max_length=max(cnv["length"],asm["length"])
    length_diff=abs(cnv["length"]-asm["length"])
    if max_length>0 and length_diff/float(max_length) > max_length_diff_frac:
        #stats_adj["badalign"][cnv["svtype"]]+=1
        #cnv["reason"]="length_diff"
        #continue
        return False,"bad_lendiff"

    if abs(cnv["start"]-asm["pos"]) > max_start_pos_diff:
        #stats_adj["badalign"][cnv["svtype"]]+=1
        #cnv["reason"]="startpos_diff"
        #continue
        return False,"bad_startdiff"

    cnv_seq=cnv["cnv_seq"]
    asm_seq=asm["event_seq"]

    alignment=pairwise2.align.globalxx(Seq(cnv_seq), Seq(asm_seq),one_alignment_only=True)[0]
    score_frac=alignment.score/max(len(cnv_seq),len(asm_seq))
    if score_frac < min_score_frac:
        #stats_adj["badalign"][cnv["svtype"]]+=1
        #cnv["reason"]="alignment_score"
        return False,"bad_alnscore"

    return True,{"aln_score":alignment.score}

#STEP1: Match CNV calls to adjacent events from asm-asm alignment

#Load putative duplication/deletion events from assembly-assembly alignment variants file into an IntervalTree
trees=collections.defaultdict(lambda: intervaltree.IntervalTree())
for var in readvcf(asm_asm_vcf_file):
    length_diff=len(var["alt"])-len(var["ref"])

    if abs(length_diff) < min_cnv_length:
        continue

    putative_svtype="DUP" if length_diff>0 else "DEL"
    stats_adj["putative"][putative_svtype]+=1

    var["putative_svtype"]=putative_svtype
    var["length"]=length_diff

    if putative_svtype=="DUP":
        var["event_seq"]=var["alt"]
    elif putative_svtype=="DEL":
        var["event_seq"]=var["ref"] #Untested
    else:
        raise Exception("Unsupported SVTYPE")

    trees[var["chr"]][(var["pos"]-padding_bp):(var["pos"]+padding_bp)]=var

cnvs=list(readvcf(cnv_vcf_file))

#Initialize CNV info and extract CNV event sequences
print("extract cnv event sequences...")
for cnv in cnvs:
    info=vcfinfo2dict(cnv["info_str"])
    start=cnv["pos"]
    end=int(info["END"])

    cnv["status"]="not_included"
    cnv["reason"]=None

    cnv["start"]=start
    cnv["end"]=end
    cnv["length"]=end-start + 1
    cnv["svtype"]=info["SVTYPE"]
    cnv["info"]=info
    cnv["good_matches_adj"]=0
    cnv["good_matches_far"]=0

    if svtype_filter!=None and len(svtype_filter) and not info["SVTYPE"] in svtype_filter:
        continue

    if cnv["length"] > max_cnv_length:
        continue

    cnv["status"]="no_match"
    cnv["reason"]="bad_other"

    cnv_seq=extractseq(ref_file,cnv["chr"],cnv["start"],cnv["end"],header=None)
    cnv["cnv_seq"]=cnv_seq
    cnv["cnv_seq_padded"]="" #extractseq(ref_file,cnv["chr"],cnv["start"]-30000,cnv["end"]+30000,header=None)

print("match cnv events (adjacent)...")
#Match each CNV
for cnv in cnvs:
    if cnv["status"]=="not_included":
        continue

    info=cnv["info"]
    start=cnv["pos"]
    end=int(info["END"])

    stats_adj["cnv_call"][info["SVTYPE"]]+=1

    overlap=trees[cnv["chr"]][start:end]

    if len(overlap)==0:
        cnv["reason"]="bad_no_overlap"
        stats_adj["bad_no_overlap"][cnv["svtype"]]+=1
        continue

    """if len(overlap)>1:
        #TODO: Implement handling of multiple hits
        cnv["reason"]="multihit"
        stats_adj["multihit"][cnv["svtype"]]+=1
        continue"""

    #for i,asm_var in enumerate(list(overlap)):
    #asmr=list(overlap)[0].data



    asm_list=[o.data for o in overlap]

    #Check each overlapping event, starting with the ones with the closest position to the CNV call start

    no_match_reason=None
    has_matched=False

    for index, asm in enumerate(sorted(asm_list,key=lambda a: abs(a["pos"]-cnv["pos"]))):
        match_ok,info=evalmatch(cnv,asm)
        if no_match_reason==None: #When no match is found, reason for being unmatched is set to the failure reason for the closest asm-asm event
            no_match_reason=info
        print(index,abs(asm["pos"]-cnv["pos"]),"cnv[",cnv["chr"],cnv["start"],cnv["length"],"] vs asm[",asm["chr"],asm["pos"],asm["length"],"]")

        if match_ok:
            cnv["good_matches_adj"]+=1

        if match_ok and not has_matched:
            cnv["matched_asm_var"]=asm
            cnv["aln_score"]=info["aln_score"]
            stats_adj["matched"][cnv["svtype"]]+=1

            cnv["reason"]="none"

            if cnv["start"]==asm["pos"] and cnv["cnv_seq"]==asm["event_seq"]:
                stats_adj["perfect"][cnv["svtype"]]+=1
                cnv["status"]="match_adj_perfect"
            else:
                cnv["status"]="match_adj"

            #if len(asm_list)>1:
            #    cnv["reason"]="multihit-%d"%

            has_matched=True

    if len(asm_list)>1:
        stats_adj["multi_overlap"][cnv["svtype"]]+=1


    if cnv["status"]=="no_match":
        cnv["reason"]=no_match_reason

        stats_adj[no_match_reason][cnv["svtype"]]+=1


    #for aln in alignments:
    #    print(aln)
        #print(pairwise2.format_alignment(*aln))

    #Implemented for duplications only for now..
    #assert(putative_svtype==var["svtype"])


print("writing results...")

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

h=open("%s/cnvlink_adjacent.fasta"%output_dir,"w")
for var in cnvs:
    if var["status"]=="adjacent_perfect" or var["status"]=="adjacent":
        h.write(">%s\n%s\n"%(var["id"],var["cnv_seq_padded"]))
h.close()

h=open("%s/cnvlink_adjacent_out.tsv"%output_dir,"w")
h.write("#generated with cnvlink %s\n"%version)
h.write("#cmd: %s\n"%" ".join(sys.argv))
h.write("#cnv_id\tstatus\tfail_reason\ttype\tgood_matches_adj\tgood_matches_far\taln_score\tcnv_chr\tcnv_start\tcnv_length\tasm_chr\tasm_start\tasm_length\tcnv_seq\tasm_seq\n")
for var in sorted(cnvs,key=lambda v: v["aln_score"] if "aln_score" in v else -1,reverse=True):
    if var["status"]=="not_included" or var["status"]=="no_match":
        continue
    h.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n"%(
                                                             var["id"],
                                                             var["status"],
                                                             var["reason"],
                                                             var["svtype"],
                                                             var["good_matches_adj"],
                                                             var["good_matches_far"],
                                                             var["aln_score"] if "aln_score" in var else -1,
                                                             var["chr"],
                                                             var["start"],
                                                             var["length"],
                                                             var["matched_asm_var"]["chr"] if "matched_asm_var" in var else "*",
                                                             var["matched_asm_var"]["pos"] if "matched_asm_var" in var else -1,
                                                             var["matched_asm_var"]["length"] if "matched_asm_var" in var else -1,
                                                             var["cnv_seq"] if "cnv_seq" in var else "*",
                                                             var["matched_asm_var"]["event_seq"] if "matched_asm_var" in var else "*"))
h.close()




print(stats_adj)
h=open("%s/cnvlink_stats_adj.json"%output_dir,"w")
h.write(json.dumps(stats_adj,sort_keys=True))
h.close()
