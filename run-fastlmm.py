from fastlmm.association import single_snp
bed_fn = "../filtered.lmm"
pheno_fn = "warpedlmm_pheno.txt"
cov_fn = "../warpedlmm.covars"
output_name = "fastlmm.out.assoc"
results_df = single_snp(bed_fn,  pheno_fn, covar=cov_fn, leave_out_one_chrom=0, output_file_name = output_name)
