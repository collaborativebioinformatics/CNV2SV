# CNVLINK
# CNV2SV project, BCM SV Hackathon October 2020
# M.S.- Mon, October 12 2020

#Prerequisites
#Python >= 3 (Tested on 3.8)
#conda install intervaltree
#conda install mappy
#conda install pyfaidx

#TODO: Confirm results
#TODO: Document output file fields (.tsv)
#TODO: Create conda package, set dependency versions

#CAUTION: How to deal with repeats called as CNVs
#CAUTION: Deletions analysis not implemented yet
#CAUTION: Output in the all_good_matches_adj and all_good_matches_far currently has the format [chr]:[start]-[length]

import os
import sys
import intervaltree
import collections
import json
import mappy
import pyfaidx

version="1.0"

if len(sys.argv)<5:
    print("cnvlink %s"%version)
    print("Usage: cnvlink.py cnv.vcf sv.vcf reference.fa output_directory_name")
    exit()

cnv_vcf_file=sys.argv[1]       #CNV calls output vcf file (calls made on GRCH38 with chm13 reads)
asm_asm_vcf_file=sys.argv[2]   #dipcall output vcf file (CHM13 vs GRCH38 assembly-assembly alignment)
ref_file=sys.argv[3]           #Reference file (GRCH38)
output_dir=sys.argv[4].strip() #Output directory to store the resulting files

svtype_filter=["DUP"]

#Parameters for confirming CNVs where duplications / deletions are adjacent
padding_bp=1000           #Padding around CNV calls added to starts/end position when recorded in intervaltree for detecting adjacent matches, larger value allows for larger distance between called CNV and asm-asm event
min_cnv_length=50         #Minimum bp for a variant called in asm_asm_vcf_file to be considered as a potential duplication/deletion event
max_cnv_length=1000*250   #Maximum length of CNVs to be included in the analysis (performance)
#max_length_diff_frac=0.2 #OBSOLETE/DROPPED Maximum length difference when matching called CNV to asm-asm events (in terms of fraction of length of longer sequence)
max_start_pos_diff=1000   #Maximum difference in starting position for adjacent duplication events
min_score_frac=0.8        #Minimum matching characters of CNV sequence within asm-asm event sequence

stats_adj={
       "putative":      {"DUP":0,"DEL":0},  #Total number of asm-asm putative dup/del events
       "cnv_call":      {"DUP":0,"DEL":0},  #Total number of CNV calls
       "multi_overlap": {"DUP":0,"DEL":0},  #CNVs that overlap multiple asm-asm events (currently, these CNVs are skipped!)
       "matched_adj":   {"DUP":0,"DEL":0},  #CNVs with a corresponding adjacent asm-asm event that fulfills all criteria
       "perfect_adj":   {"DUP":0,"DEL":0},  #CNVs with a corresponding adjacent asm-asm event that fulfills all criteria and exactly match in terms of start position and length
       "matched_far":   {"DUP":0,"DEL":0},  #CNVs with a corresponding non-adjacent asm-asm event that fulfills all criteria
       "perfect_far":   {"DUP":0,"DEL":0},  #CNVs with a corresponding non-adjacent asm-asm event that fulfills all criteria and the alignment match length is equal to all bases from the CNV
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

def evalmatch(cnv,asm,ignore_start_pos_diff=False):
    #Matched e.g. a CNV duplication call with a deletion event
    if asm["putative_svtype"] != cnv["svtype"]:
        #cnv["reason"]="nonsense"
        #stats_adj["nonsense"][cnv["svtype"]]+=1
        return False,"bad_nonsense"

    """#Length difference of reported event too long
    max_length=max(cnv["length"],asm["length"])
    length_diff=abs(cnv["length"]-asm["length"])
    if max_length>0 and length_diff/float(max_length) > max_length_diff_frac:
        #stats_adj["badalign"][cnv["svtype"]]+=1
        #cnv["reason"]="length_diff"
        #continue
        return False,"bad_lendiff"""

    if not ignore_start_pos_diff:
        if abs(cnv["start"]-asm["pos"]) > max_start_pos_diff:
            #stats_adj["badalign"][cnv["svtype"]]+=1
            #cnv["reason"]="startpos_diff"
            #continue
            return False,"bad_startdiff"

    cnv_seq=cnv["cnv_seq"]
    asm_seq=asm["event_seq"]

    aligner=mappy.Aligner(seq=asm_seq)
    alignments=list(aligner.map(cnv_seq))
    best_score=0
    if len(alignments)==0:
        return False,"bad_alnscore"

    best_score=max([aln.mlen for aln in alignments])
    best_aln=[aln for aln in alignments if aln.mlen==best_score][0]

    score_frac=best_score/float(len(cnv_seq))
    if score_frac < min_score_frac:
        #stats_adj["badalign"][cnv["svtype"]]+=1
        #cnv["reason"]="alignment_score"
        return False,"bad_alnscore"

    return True,{"alignment":best_aln}

#STEP1: Match CNV calls to adjacent events from asm-asm alignment

#Load putative duplication/deletion events from assembly-assembly alignment variants file into an IntervalTree
all_asm=list(readvcf(asm_asm_vcf_file))
trees=collections.defaultdict(lambda: intervaltree.IntervalTree())
index=0
for asm in all_asm:
    length_diff=len(asm["alt"])-len(asm["ref"])

    asm["included"]=False
    if abs(length_diff) < min_cnv_length:
        continue

    putative_svtype="DUP" if length_diff>0 else "DEL"
    stats_adj["putative"][putative_svtype]+=1

    if svtype_filter!=None and len(svtype_filter) and not putative_svtype in svtype_filter:
        continue

    asm["included"]=True
    asm["putative_svtype"]=putative_svtype
    asm["length"]=length_diff
    asm["internal_id"]="ASM_%08d"%index
    index+=1

    if putative_svtype=="DUP":
        asm["event_seq"]=asm["alt"]
    elif putative_svtype=="DEL":
        asm["event_seq"]=asm["ref"] #Untested
    else:
        raise Exception("Unsupported SVTYPE")

    trees[asm["chr"]][(asm["pos"]-padding_bp):(asm["pos"]+padding_bp)]=asm

"""print("load cnv sequences...")
h=open("cnvs.tmp","r")
cnvs=json.loads(h.read())
h.close()"""


cnvs=list(readvcf(cnv_vcf_file))

#Initialize CNV info and extract CNV event sequences
print("extract cnv event sequences...")
ref=pyfaidx.Fasta(ref_file)

index=0
last_pct=None
for cnv in cnvs:
    info=vcfinfo2dict(cnv["info_str"])
    start=cnv["pos"]
    end=int(info["END"])

    cnv["status"]="not_included"
    cnv["reason"]=None

    cnv["start"]=start
    cnv["end"]=end
    cnv["length"]=end-start + 1
    cnv["internal_id"]="CNV_%08d"%index
    index+=1
    cnv["svtype"]=info["SVTYPE"]
    cnv["info"]=info
    cnv["good_matches_adj"]=[]
    cnv["good_matches_far"]=[]

    pct_done=int(100*index/float(len(cnvs)))
    if pct_done%5==0 and pct_done!=last_pct:
        print("%d%%"%pct_done)
        last_pct=pct_done

    if svtype_filter!=None and len(svtype_filter) and not info["SVTYPE"] in svtype_filter:
        cnv["reason"]="unsupported_svtype"
        continue

    if cnv["length"] > max_cnv_length:
        cnv["reason"]="too_long"
        continue

    cnv["status"]="no_match"
    cnv["reason"]="bad_other"
    cnv["cnv_seq"]=str(ref[cnv["chr"]][cnv["start"]-1:cnv["end"]])
    #cnv["cnv_seq_padded"]="" #extractseq(ref_file,cnv["chr"],cnv["start"]-30000,cnv["end"]+30000,header=None)


"""
cnv_temp_out=open("cnvs.tmp","w")
cnv_temp_out.write(json.dumps(cnvs))
cnv_temp_out.close()"""

print("match cnv events (adjacent)...")
#Search for adjacent matching events for each each CNV
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

    asm_overlap_list=[o.data for o in overlap]

    no_match_reason=None
    has_matched=False

    for index, asm in enumerate(sorted(asm_overlap_list,key=lambda a: abs(a["pos"]-cnv["pos"]))):
        match_ok,info=evalmatch(cnv,asm)
        if no_match_reason==None: #When no match is found, reason for being unmatched is set to the failure reason for the closest asm-asm event
            no_match_reason=info
        print(index,abs(asm["pos"]-cnv["pos"]),"cnv[",cnv["chr"],cnv["start"],cnv["length"],"] vs asm[",asm["chr"],asm["pos"],asm["length"],"]")

        if match_ok:
            cnv["good_matches_adj"].append(asm)

        if match_ok and not has_matched:
            cnv["matched_asm_var"]=asm
            cnv["alignment"]=info["alignment"]
            stats_adj["matched_adj"][cnv["svtype"]]+=1

            cnv["reason"]="none"

            if cnv["start"]==asm["pos"] and cnv["cnv_seq"]==asm["event_seq"]:
                stats_adj["perfect_adj"][cnv["svtype"]]+=1
                cnv["status"]="match_adj_perfect"
            else:
                cnv["status"]="match_adj"

            #if len(asm_list)>1:
            #    cnv["reason"]="multihit-%d"%

            has_matched=True

    if len(asm_overlap_list)>1:
        stats_adj["multi_overlap"][cnv["svtype"]]+=1


    if cnv["status"]=="no_match":
        cnv["reason"]=no_match_reason

        stats_adj[no_match_reason][cnv["svtype"]]+=1


    #for aln in alignments:
    #    print(aln)
        #print(pairwise2.format_alignment(*aln))

    #Implemented for duplications only for now..
    #assert(putative_svtype==var["svtype"])


if not os.path.exists(output_dir):
    os.mkdir(output_dir)

asm_seqs=""
asm_id_to_asm={}
for asm in all_asm:
    if not asm["included"]:
        continue
    asm_seqs+=">%s\n%s\n"%(asm["internal_id"],asm["event_seq"])
    asm_id_to_asm[asm["internal_id"]]=asm

print("match cnv events (far)...")

print("init mappy aligner...")
asm_seq_file=output_dir+"/asm.tmp.fa"
asm_seq_handle=open(asm_seq_file,"w")
asm_seq_handle.write(asm_seqs)
asm_seq_handle.close()
asm_aligner=mappy.Aligner(asm_seq_file)


# STEP 2: Look for matching non-adjacent events
print("searching for non-adjacent links...")
done_count=0
for cnv in cnvs:
    done_count+=1
    if cnv["status"]=="not_included":
        continue

    if done_count%1000==0:
        print("searching... %d/%d (%.2f%%)"%(done_count,len(cnvs),100.0*done_count/float(len(cnvs))))

    #print(cnv,len(cnv["cnv_seq"]),cnv["length"])
    alignments=sorted(asm_aligner.map(cnv["cnv_seq"]),key=lambda a: a.mlen,reverse=True)
    for i,aln in enumerate(alignments):
        #print(aln)
        #import code
        #code.interact(local=locals())

        #print(aln.ctg,asm_id_to_asm[aln.ctg]["length"],len(asm_id_to_asm[aln.ctg]["event_seq"]))

        asm=asm_id_to_asm[aln.ctg] #aln.ctg: Reference sequence name
        match_ok,info=True,""

        cnv_seq=cnv["cnv_seq"]
        asm_seq=asm["event_seq"]

        score_frac=aln.mlen/float(max(len(cnv_seq),len(asm_seq)))
        if score_frac < min_score_frac:
            match_ok,info=False,"bad_alnscore"

        #Matched e.g. a CNV duplication call with a deletion event
        if asm["putative_svtype"] != cnv["svtype"]:
            match_ok,info=False,"bad_nonsense"

        #Length difference of reported event too long
        """max_length=max(cnv["length"],asm["length"])
        length_diff=abs(cnv["length"]-asm["length"])
        if max_length>0 and length_diff/float(max_length) > max_length_diff_frac:
            match_ok,info=False,"bad_lendiff"
        """

        if match_ok:
            if not asm in cnv["good_matches_adj"]:
                cnv["good_matches_far"].append(asm)
            if cnv["status"]=="no_match":
                if len(cnv["cnv_seq"])==aln.mlen:
                    stats_adj["perfect_far"][cnv["svtype"]]+=1
                    cnv["status"]="match_far_perfect"
                else:
                    cnv["status"]="match_far"

                cnv["matched_asm_var"]=asm
                cnv["alignment"]=aln
                stats_adj["matched_far"][cnv["svtype"]]+=1

                cnv["reason"]="none"


        #else:
        #    print("BAD MATCH",info,aln.mlen,len(cnv_seq),len(asm_seq),aln.mlen/float(max(len(cnv_seq),len(asm_seq))))

    """for asm in all_asm:
        done_comparisons+=1
        if not asm["included"] or asm["putative_svtype"]!=cnv["svtype"]:
            continue
        match_ok,info=evalmatch(cnv,asm,ignore_start_pos_diff=True)
        if match_ok:
            print(abs(asm["pos"]-cnv["pos"]),"cnv[",cnv["chr"],cnv["start"],cnv["length"],"] vs asm[",asm["chr"],asm["pos"],asm["length"],"]")

        if done_comparisons%100000==0:
            print("FAR %d/%d (%.2f%%)"%(done_comparisons,total_comparisons,100*done_comparisons/float(total_comparisons)))"""





print("writing results...")


"""
h=open("%s/cnvlink_adjacent.fasta"%output_dir,"w")
for var in cnvs:
    if var["status"]=="adjacent_perfect" or var["status"]=="adjacent":
        h.write(">%s\n%s\n"%(var["id"],var["cnv_seq_padded"]))
h.close()"""

def fmt_all_matches(matches):
    return ",".join(["%s:%d-%d"%(m["chr"],m["pos"],m["length"]) for m in matches])

def write_results_tsv(filename,filter_bad=True):
    h=open(filename,"w")
    h.write("#generated with cnvlink %s\n"%version)
    h.write("#cmd: %s\n"%" ".join(sys.argv))
    h.write("#cnv_id\tstatus\tfail_reason\ttype\tgood_matches_adj\tgood_matches_far\talignment_matches\tcnv_chr\tcnv_start\tcnv_length\tasm_chr\tasm_start\tasm_length\tall_good_matches_adj\tall_good_matches_far\taln_cigar\taln_length\taln_NM\taln_mapq\tcnv_seq\tasm_seq\n")
    for var in sorted(cnvs,key=lambda v: v["alignment"].mlen if "alignment" in v else -1,reverse=True):
        if filter_bad and (var["status"]=="not_included" or var["status"]=="no_match"):
            continue
        h.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n"%(
                                                                 var["id"],
                                                                 var["status"],
                                                                 var["reason"],
                                                                 var["svtype"],
                                                                 len(var["good_matches_adj"]),
                                                                 len(var["good_matches_far"]),
                                                                 var["alignment"].mlen if "alignment" in var else -1,
                                                                 var["chr"],
                                                                 var["start"],
                                                                 var["length"],
                                                                 var["matched_asm_var"]["chr"] if "matched_asm_var" in var else "*",
                                                                 var["matched_asm_var"]["pos"] if "matched_asm_var" in var else -1,
                                                                 var["matched_asm_var"]["length"] if "matched_asm_var" in var else -1,
                                                                 fmt_all_matches(var["good_matches_adj"]),
                                                                 fmt_all_matches(var["good_matches_far"]),
                                                                 var["alignment"].cigar_str if "alignment" in var else "",
                                                                 var["alignment"].blen if "alignment" in var else -1,
                                                                 var["alignment"].NM if "alignment" in var else -1,
                                                                 var["alignment"].mapq if "alignment" in var else -1,
                                                                 var["cnv_seq"] if "cnv_seq" in var else "*",
                                                                 var["matched_asm_var"]["event_seq"] if "matched_asm_var" in var else "*"))
    h.close()

write_results_tsv("%s/cnvlink_out.tsv"%output_dir,filter_bad=True)
write_results_tsv("%s/cnvlink_out_all.tsv"%output_dir,filter_bad=False)


print(stats_adj)
h=open("%s/cnvlink_stats.json"%output_dir,"w")
h.write(json.dumps(stats_adj,sort_keys=True))
h.close()
