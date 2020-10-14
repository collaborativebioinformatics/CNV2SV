# CNVLINK
# CNV2SV project, BCM SV Hackathon October 2020
# M.S.- Mon, October 12 2020

import sys

from_bed=open("from.bed","w")
from_bed.write("chr\tstart\tend\tvalue1\n")

to_bed=open("to.bed","w")
to_bed.write("chr\tstart\tend\tvalue1\n")

best_match_only=True

for line_ in open(sys.argv[1],"r").readlines()[3:]:
    line=line_.strip()
    cnv_id,status,fail_reason,type,good_matches_adj,good_matches_far,alignment_matches,cnv_chr,cnv_start,cnv_length,asm_chr,asm_start,asm_length,all_good_matches_adj,all_good_matches_far,aln_cigar,aln_length,aln_NM,aln_mapq,cnv_seq,asm_seq=line.split("\t")
    cnv_start=int(cnv_start)
    cnv_length=int(cnv_length)
    if status=="match_far" or status=="match_far_perfect":
        matches=all_good_matches_far.split(",")
        if best_match_only:
            matches=[matches[0]]
        for match in matches:
            chr_from=cnv_chr
            chr_to=match.split(":")[0]
            pos_to=match.split(":")[1].split("-")[0]
            from_bed.write("%s\t%s\t%s\t1.0\n"%(cnv_chr,cnv_start-10**7.0/2.0,cnv_start+10**7.0/2.0))
            to_bed.write("%s\t%s\t%s\t1.0\n"%(chr_to,int(pos_to)-10**6.3/2.0,int(pos_to)+10**6.3/2.0))

    if status=="match_adj" or status=="match_adj_perfect":
        matches=all_good_matches_adj.split(",")
        if best_match_only:
            matches=[matches[0]]
        for match in matches:
            chr_from=cnv_chr
            chr_to=match.split(":")[0]
            pos_to=match.split(":")[1].split("-")[0]
            from_bed.write("%s\t%s\t%s\t1.0\n"%(cnv_chr,cnv_start-10**7.0/2.0,cnv_start+10**7.0/2.0))
            to_bed.write("%s\t%s\t%s\t1.0\n"%(chr_to,int(pos_to)-10**6.3/2.0,int(pos_to)+10**6.3/2.0))

from_bed.close()
to_bed.close()
