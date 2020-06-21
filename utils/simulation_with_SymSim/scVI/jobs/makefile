sim_data_script := symsim_data_transform.R
impute_script := impute.sh

generate:
	Rscript --vanilla ../${sim_data_script}

# impute:
# 	for i in 1..20; do\
#      for p in umi nonumi; do\
#          python ../${impute_script} --seed $${i} --platform $${p} >$${i}_$${p}.log 2>&1 &; \
#       done; \
#   done

impute:
	./${impute_script}

rm_symsim:
	rm ../data/symsim/*

rm_impute:
	rm ../../simutool/jobs/scvi_result/*

rm_log:
	rm ./*.log
