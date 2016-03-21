args = commandArgs(T)
purity_summary_table = toString(args[1])

TEMPLATE = "template_rho_and_psi.txt"

dat = read.table(purity_summary_table, header=T, stringsAsFactors=F)
rho_psi_template = read.table(TEMPLATE, header=T, stringsAsFactors=F)

for (i in 1:nrow(dat)) {
	samplename = dat[i, "tumor"]
	rho_psi_template["FRAC_GENOME", "rho"] = dat[i, "purity_consensus"]
	write.table(rho_psi_template, file=paste0(samplename, "_rho_and_psi.txt"), row.names=T, quote=F, sep="\t")
}
