using MAT
using DelimitedFiles

fname = "yfin_WT_NaGain_1Hz"

writedlm("$(fname).tsv", matread("$(fname).mat")["yfinal"], '\t')
