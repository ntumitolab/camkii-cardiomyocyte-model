using MAT
using DelimitedFiles
using CaMKIIModel

fname = "yfin_WT_NaGain_1Hz"

writedlm("$(fname).tsv", matread("$(fname).mat")["yfinal"], '\t')
