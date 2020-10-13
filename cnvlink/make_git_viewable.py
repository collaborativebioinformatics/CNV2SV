

def make_git_viewable(dir):
    tsv=open(dir+"/cnvlink_out.tsv","r").readlines()
    tsv=tsv[2:]
    tsv_git=open(dir+"/cnvlink_out_git.tsv","w")
    for line_ in tsv:
        line=line_.strip()
        parts=line.split("\t")
        assert(len(parts)==21)
        tsv_git.write("\t".join(parts[:19])+"\n")

make_git_viewable("cnvlink_out_freec")
make_git_viewable("cnvlink_out_parliament_cnvnator")
make_git_viewable("cnvlink_out_parliament_combined")
make_git_viewable("cnvlink_out_parliament_subset_cnvnator")
