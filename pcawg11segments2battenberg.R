
args = commandArgs(T)
infile = toString(args[1])
outfile = toString(args[2])

new = read.table(infile, header=T, stringsAsFactors=F)
out = cbind(new[, c("chromosome", "start", "end")], matrix(NA, ncol=4, nrow=nrow(new)), new[, c("major_cn", "minor_cn")], rep(1, nrow(new)), matrix(NA, ncol=(7 + 10*5), nrow=nrow(new)))

  colnames(out) = c("chr","startpos","endpos","BAF","pval","LogR","ntot",
                            "nMaj1_A","nMin1_A","frac1_A","nMaj2_A","nMin2_A","frac2_A","SDfrac_A","SDfrac_A_BS","frac1_A_0.025","frac1_A_0.975",
                            "nMaj1_B","nMin1_B","frac1_B","nMaj2_B","nMin2_B","frac2_B","SDfrac_B","SDfrac_B_BS","frac1_B_0.025","frac1_B_0.975",
                            "nMaj1_C","nMin1_C","frac1_C","nMaj2_C","nMin2_C","frac2_C","SDfrac_C","SDfrac_C_BS","frac1_C_0.025","frac1_C_0.975",
                            "nMaj1_D","nMin1_D","frac1_D","nMaj2_D","nMin2_D","frac2_D","SDfrac_D","SDfrac_D_BS","frac1_D_0.025","frac1_D_0.975",
                            "nMaj1_E","nMin1_E","frac1_E","nMaj2_E","nMin2_E","frac2_E","SDfrac_E","SDfrac_E_BS","frac1_E_0.025","frac1_E_0.975",
                            "nMaj1_F","nMin1_F","frac1_F","nMaj2_F","nMin2_F","frac2_F","SDfrac_F","SDfrac_F_BS","frac1_F_0.025","frac1_F_0.975")
write.table(out, file=outfile, row.names=F, quote=F, sep="\t")
