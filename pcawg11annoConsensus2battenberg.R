args = commandArgs(T)
infile_segments = args[1]
infile_purity = args[2]
outfile_segments = args[3]
outfile_purity = args[4]

TEMPLATE = "~/repo/icgcConsensusCN_to_BB/template_rho_and_psi.txt"

set_colnames = function(df) {
	colnames(df) = c("chr","startpos","endpos","BAF","pval","LogR","ntot",
                            "nMaj1_A","nMin1_A","frac1_A","nMaj2_A","nMin2_A","frac2_A","SDfrac_A","SDfrac_A_BS","frac1_A_0.025","frac1_A_0.975",
                            "nMaj1_B","nMin1_B","frac1_B","nMaj2_B","nMin2_B","frac2_B","SDfrac_B","SDfrac_B_BS","frac1_B_0.025","frac1_B_0.975",
                            "nMaj1_C","nMin1_C","frac1_C","nMaj2_C","nMin2_C","frac2_C","SDfrac_C","SDfrac_C_BS","frac1_C_0.025","frac1_C_0.975",
                            "nMaj1_D","nMin1_D","frac1_D","nMaj2_D","nMin2_D","frac2_D","SDfrac_D","SDfrac_D_BS","frac1_D_0.025","frac1_D_0.975",
                            "nMaj1_E","nMin1_E","frac1_E","nMaj2_E","nMin2_E","frac2_E","SDfrac_E","SDfrac_E_BS","frac1_E_0.025","frac1_E_0.975",
                            "nMaj1_F","nMin1_F","frac1_F","nMaj2_F","nMin2_F","frac2_F","SDfrac_F","SDfrac_F_BS","frac1_F_0.025","frac1_F_0.975")
	return(df)
}

samplename = unlist(strsplit(basename(infile_segments), ".", fixed=T))[1]
dat = read.table(infile_segments, header=T, stringsAsFactors=F)
# Removing segments that should not have been released
dat = dat[dat$total_cn >=0,]

# mandatory star 3 segments first
dat_sub = dat[dat$star == 3,]
out1 = cbind(dat_sub[, c("chromosome", "start", "end")], matrix(NA, ncol=4, nrow=nrow(dat_sub)), dat_sub[, c("major_cn", "minor_cn")], rep(1, nrow(dat_sub)), matrix(NA, ncol=(7 + 10*5), nrow=nrow(dat_sub)))
out1 = set_colnames(out1)


dat_sub = dat[dat$star != 3,]

# Fill in rest with Battenberg calls - if available
if (any(!is.na(dat$battenberg_nMaj1_A))) {
	subclonal_method = "battenberg"
	cat(paste0(samplename, "\t", subclonal_method, "\n"))
	out2 = cbind(dat_sub[, c("chromosome", "start", "end")], matrix(NA, ncol=4, nrow=nrow(dat_sub)), dat_sub[, c("battenberg_nMaj1_A", "battenberg_nMin1_A", "battenberg_frac1_A", "battenberg_nMaj2_A", "battenberg_nMin2_A", "battenberg_frac2_A", "battenberg_SDfrac_A", "battenberg_SDfrac_A_BS", "battenberg_frac1_A_0.025", "battenberg_frac1_A_0.975")], matrix(NA, ncol=(10*5), nrow=nrow(dat_sub)))

} else if (any(!is.na(dat$sclust_nMaj1_A))) {
# Or alternatively Sclust - if available
	subclonal_method = "sclust"
	cat(paste0(samplename, "\t", subclonal_method, "\n"))
	out2 = cbind(dat_sub[, c("chromosome", "start", "end")], matrix(NA, ncol=4, nrow=nrow(dat_sub)), dat_sub[, c("sclust_nMaj1_A", "sclust_nMin1_A", "sclust_frac1_A", "sclust_nMaj2_A", "sclust_nMin2_A", "sclust_frac2_A")], matrix(NA, ncol=(4 + 10*5), nrow=nrow(dat_sub)))
} else {
	subclonal_method = "purity"
	cat(paste0(samplename, "\t", "consensus", "\n"))
	out2 = cbind(dat_sub[, c("chromosome", "start", "end")], matrix(NA, ncol=4, nrow=nrow(dat_sub)), dat_sub[, c("major_cn", "minor_cn")], rep(1, nrow(dat_sub)), matrix(NA, ncol=(7 + 10*5), nrow=nrow(dat_sub)))
}
out2 = set_colnames(out2)

# Combine and sort
out = as.data.frame(rbind(out1, out2))
out = out[with(out, order(chr, startpos)),]

write.table(out, file=outfile_segments, row.names=F, quote=F, sep="\t")

# Write out the purity
rho_psi_template = read.table(TEMPLATE, header=T, stringsAsFactors=F)
all_purities = read.table(infile_purity, header=T, stringsAsFactors=F)
purity = all_purities[all_purities$samplename==samplename, colnames(all_purities)==subclonal_method]

rho_psi_template["FRAC_GENOME", "rho"] = purity
write.table(rho_psi_template, file=outfile_purity, row.names=T, quote=F, sep="\t")
