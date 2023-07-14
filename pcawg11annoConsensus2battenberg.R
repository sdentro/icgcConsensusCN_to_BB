args = commandArgs(T)
infile_segments = args[1]
infile_purity = args[2]
outfile_segments = args[3]
outfile_purity = args[4]
subset_run = as.integer(args[5])

TEMPLATE = "~/repo/icgcConsensusCN_to_BB/template_rho_and_psi.txt"
LEVELS_TO_REPLACE = c("e", "f", "g", "h", "i")

samplename = unlist(strsplit(basename(infile_segments), ".", fixed=T))[1]
dat = read.table(infile_segments, header=T, stringsAsFactors=F)
# Removing segments that should not have been released
sel = dat$total_cn >=0
sel[is.na(sel)] = F
dat = dat[sel,]

# Perform subsetting depending on the run
if (subset_run==1) {
	keep_consensus = dat$level %in% c("a", "b", "c", "d")
	replace_cn = F
} else if (subset_run==2) {
	keep_consensus = dat$star %in% c(3,2)
	replace_cn = F
} else if (subset_run==3) {
	keep_consensus = dat$star %in% c(3,2)
	replace_cn = T
} else if (subset_run==4) {
	keep_consensus = rep(T, nrow(dat))
	replace_cn = F
} else if (subset_run==5) {
	keep_consensus = rep(T, nrow(dat))
	replace_cn = T
} else if (subset_run==6) {
	keep_consensus = dat$star %in% c(3)
        replace_cn = F
} else {
	print("Please supply 1-6 to denote the combination of segments to be selected")
	q(save="no")
}



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

# Subset by selection
all_purities = read.table(infile_purity, header=T, stringsAsFactors=F)
dat_sub = dat[keep_consensus,]
out1 = cbind(dat_sub[, c("chromosome", "start", "end")], matrix(NA, ncol=4, nrow=nrow(dat_sub)), dat_sub[, c("major_cn", "minor_cn")], rep(1, nrow(dat_sub)), matrix(NA, ncol=(7 + 10*5), nrow=nrow(dat_sub)))

if (replace_cn) {
	to_replace = dat_sub$level %in% LEVELS_TO_REPLACE

	# Fill in rest with Battenberg calls - if available
	if (sum(to_replace, na.rm=T) > 0 & any(!is.na(dat_sub$battenberg_nMaj1_A)) & !is.na(all_purities[all_purities$samplename==samplename, "battenberg"])) {
		subclonal_method = "battenberg"
		cat(paste0(samplename, "\t", subset_run, "\t", subclonal_method, "\n"))
		dat_replace = cbind(dat_sub[to_replace, c("chromosome", "start", "end")], 
				    matrix(NA, ncol=4, nrow=sum(to_replace, na.rm=T)), 
				    dat_sub[to_replace, c("battenberg_nMaj1_A", "battenberg_nMin1_A", "battenberg_frac1_A", "battenberg_nMaj2_A", "battenberg_nMin2_A", "battenberg_frac2_A", "battenberg_SDfrac_A", "battenberg_SDfrac_A_BS", "battenberg_frac1_A_0.025", "battenberg_frac1_A_0.975")], 
				    matrix(NA, ncol=(10*5), nrow=sum(to_replace, na.rm=T)))
		dat_replace = set_colnames(dat_replace)
		out1[to_replace,] = dat_replace

	
	} else if (any(!is.na(dat$sclust_nMaj1_A)) & !is.na(all_purities[all_purities$samplename==samplename, "sclust"])) {
	# Or alternatively Sclust - if available
		subclonal_method = "sclust"
		cat(paste0(samplename, "\t", subset_run, "\t", subclonal_method, "\n"))
		dat_replace = cbind(dat_sub[to_replace, c("chromosome", "start", "end")], 
				    matrix(NA, ncol=4, nrow=sum(to_replace, na.rm=T)), 
				    dat_sub[to_replace, c("sclust_nMaj1_A", "sclust_nMin1_A", "sclust_frac1_A", "sclust_nMaj2_A", "sclust_nMin2_A", "sclust_frac2_A")], 
				    matrix(NA, ncol=(4 + 10*5), nrow=sum(to_replace, na.rm=T)))
		dat_replace = set_colnames(dat_replace)
		out1[to_replace,] = dat_replace

	} else if (length(to_replace) > 0) {
		subclonal_method = "purity"
		cat(paste0(samplename, "\t", subset_run, "\t", "not_replaced - consensus", "\n"))
	} else {
		subclonal_method = "purity"
		cat(paste0(samplename, "\t", subset_run, "\t", "nothing_to_replace - consensus", "\n"))
	}
} else {
	# Nothing to replace, use consensus regardless
	subclonal_method = "purity"
	cat(paste0(samplename, "\t", subset_run, "\t", "consensus", "\n"))
}
out1 = set_colnames(out1)

# Combine and sort
out = out1[with(out1, order(chr, startpos)),]
write.table(out, file=outfile_segments, row.names=F, quote=F, sep="\t")

# Write out the purity
rho_psi_template = read.table(TEMPLATE, header=T, stringsAsFactors=F)
all_purities = read.table(infile_purity, header=T, stringsAsFactors=F)
purity = all_purities[all_purities$samplename==samplename, colnames(all_purities)==subclonal_method]
rho_psi_template["FRAC_GENOME", "rho"] = purity
write.table(rho_psi_template, file=outfile_purity, row.names=T, quote=F, sep="\t")
