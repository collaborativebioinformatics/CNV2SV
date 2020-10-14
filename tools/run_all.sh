python cnvlink.py input_cnv_calls/parliament2_subset/chm13_markdup3.cnvnator.vcf ../../data/chm13_grch38.pair.vcf ../../liftover/GRCh38.no_alt_analysis_set.fa cnvlink_out_parliament_subset_cnvnator
python cnvlink.py input_cnv_calls/parliament2/PCRfree.cnvnator.vcf ../../data/chm13_grch38.pair.vcf ../../liftover/GRCh38.no_alt_analysis_set.fa cnvlink_out_parliament_cnvnator
python cnvlink.py input_cnv_calls/parliament2/PCRfree.combined.genotyped.vcf ../../data/chm13_grch38.pair.vcf ../../liftover/GRCh38.no_alt_analysis_set.fa cnvlink_out_parliament_combined
python cnvlink.py input_cnv_calls/FREEC_out/chm13_markdup3.bam_CNVs.vcf ../../data/chm13_grch38.pair.vcf ../../liftover/GRCh38.no_alt_analysis_set.fa cnvlink_out_freec
