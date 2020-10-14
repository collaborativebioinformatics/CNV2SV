import sys
#
if len(sys.argv)<2:
    print("Converts FREEC output to .vcf file readable by cnvlink.py")
    print("usage example: python freec_bed_to_vcf.py FREEC_out/chm13_markdup3.bam_CNVs > FREEC_out/chm13_markdup3.bam_CNVs.vcf")
    exit()

id=0
for line in open(sys.argv[1]):
    chr_num,start,end,cn,event_freec=line.strip().split("\t")

    if event_freec=="gain":
        event="DUP"
    elif event_freec=="loss":
        event="DEL"
    else:
        raise

    info="END=%s;SVTYPE=%s;SVLEN=%s;SVMETHOD=FREEC"%(end,event,int(end)-int(start))

    chr="chr"+chr_num
    print("%s\t%s\tFREEC_%05d\tN\t<%s>\t*\t*\t%s\t*\t*"%(chr,start,id,event,info))
    id+=1
