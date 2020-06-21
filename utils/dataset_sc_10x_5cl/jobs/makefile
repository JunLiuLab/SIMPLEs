rscript := impute_cluster_de.R
jobscript := 10xGenomics.sh
impt_script := impute.sh
clust_script := clustering.sh
no_impt_jobscript := noimpt_10xGenomics.sh
de_after_impt_jobscript := de_after_impt_10x.sh


other_methods := control scrabble viper saver rmagic scimpute
de_other_methods := control scrabble viper saver rmagic scimpute scvi
method ?= control
tests := t wilcox
thresholds := 1.0 1.2 1.4 1.5 2.0 2.1 2.2 2.4 2.8 3.0

mydir:
	-mkdir ../result
	-mkdir ../impt

test:
	Rscript --vanilla ../${rscript} -m ${method} -t wilcox

tsub:
	sbatch ${jobscript} ${method} 2 1 2 t

other_methods:
	for m in ${other_methods}; do \
      for t in ${tests}; do \
          sbatch ${jobscript} $${m} 2 1 $${t} ;\
      done; \
  done

tsimples:
	for k in 10; do \
     for m in 1; do \
        for t in ${tests}; do \
            sbatch ${jobscript} simple $${k} $${m} $${t};\
        done; \
     done; \
  done


simples:
	for k in 10 15 20; do \
     for m in 1 4 10; do \
        for t in ${tests}; do \
            sbatch ${jobscript} simple $${k} $${m} $${t};\
        done; \
     done; \
  done

b2simples:
	for k in 15 20; do \
     for m in 4 10; do \
        for t in ${tests}; do \
            sbatch ${jobscript} simple $${k} $${m} $${t};\
        done; \
     done; \
  done

b1simples:
	for k in 10; do \
     for m in 10; do \
        for t in ${tests}; do \
            sbatch ${jobscript} simple $${k} $${m} $${t};\
        done; \
     done; \
  done


local_scvi:
	for t in ${tests}; do \
	    Rscript --vanilla ../${rscript} -m scvi -t $${t}; \
  done

sub_scvi:
	for t in ${tests}; do \
	    sbatch ${no_impt_jobscript} scvi 2 1 $${t}; \
  done

test_de_simples:
	for k in 10 ; do \
     for m in 1 ; do \
        for t in t; do \
            for thres in 1.2; do \
                sbatch ${de_after_impt_jobscript} simple $${k} $${m} $${t} $${thres};\
            done; \
        done; \
     done; \
  done


de_simples:
	for k in 10 15 20; do \
     for m in 1 4 10; do \
        for t in ${tests}; do \
            for thres in ${thresholds}; do \
                sbatch ${de_after_impt_jobscript} simple $${k} $${m} $${t} $${thres};\
            done; \
        done; \
     done; \
  done

test_de_other_methods:
	sbatch ${de_after_impt_jobscript} ${method} 2 1 t 2.0;\


de_other_methods:
	for m in ${other_methods}; do \
      for t in ${tests}; do \
					for thres in ${thresholds}; do \
							sbatch ${de_after_impt_jobscript} $${m} 2 1 $${t} $${thres};\
					done; \
      done; \
  done

de_single_other:
	for t in ${tests}; do \
			for thres in ${thresholds}; do \
					sbatch ${de_after_impt_jobscript} ${method} 2 1 $${t} $${thres};\
			done; \
  done

impt_other_methods:
	sbatch ${impt_script} ${method} 2 1 wilcox

clust_single_method:
	sbatch ${clust_script} ${method} 2 1 wilcox


clean:
	rm *.err
	rm *.out
