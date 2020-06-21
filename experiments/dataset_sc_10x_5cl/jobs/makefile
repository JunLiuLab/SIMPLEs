rscript := impute_cluster_de.R

other_methods := control scrabble viper saver rmagic scimpute
method ?= control
tests := t wilcox
thresholds := 1.0 1.2 1.4 1.5 2.0 2.1 2.2 2.4 2.8 3.0

mydir:
	-mkdir ../result
	-mkdir ../impt

test:
	Rscript --vanilla ../${rscript} -m ${method} -t wilcox

other_methods:
	for m in ${other_methods}; do \
      for t in ${tests}; do \
          sbatch ${rscript} $${m} 2 1 $${t} ;\
      done; \
  done

simples:
	for k in 10 15 20; do \
     for m in 1 4 10; do \
        for t in ${tests}; do \
            sbatch ${rscript} simple $${k} $${m} $${t};\
        done; \
     done; \
  done

# This needs to run the de_after_impute function in ../impute_cluster_de.R.
local_scvi:
	for t in ${tests}; do \
	    Rscript --vanilla ../${rscript} -m scvi -t $${t}; \
  done

clean:
	rm *.err
	rm *.out
