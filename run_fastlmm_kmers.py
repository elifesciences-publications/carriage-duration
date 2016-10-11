from fastlmm.association import single_snp
import os

array_num = os.environ['LSB_JOBINDEX']

bed_fn = "kmers_" + str(array_num)
K0_fn = "../filtered.lmm"
pheno_fn = "warpedlmm_pheno.txt"
cov_fn = "../warpedlmm.covars"
cache_fn = "kmer_lmm_cache"
output_name = "fastlmm.out." + str(array_num) + ".assoc"
results_df = single_snp(bed_fn,  pheno_fn, K0=K0_fn, cache_file = cache_fn, covar=cov_fn, leave_out_one_chrom=0, output_file_name = output_name)
