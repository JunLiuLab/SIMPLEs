rscript := impute_tsne_clustering.R
jobscript := impute_tsne_clustering.sh
sjobscript := impute_tsne_clustering_simple.sh
other_timer_script := others_timer.sh
other_methods := control rmagic saver
others_of_timer := rmagic saver scimpute scrabble viper
method ?= control
cnm ?= 12000

sub_others:
	$(foreach var, ${other_methods}, sbatch ${jobscript} ${var} 2 1 ; )

sub_method:
	sbatch ${jobscript} ${method} 2 1

tsub_simples:
	for k in 20; do \
     for m in 8; do \
         sbatch ${sjobscript} simple $${k} $${m} ;\
     done; \
  done


sub_simples:
	for k in 5 10 15 20; do \
      for m in 1 2 8 12 22; do \
          sbatch ${sjobscript} simple $${k} $${m}; \
      done; \
  done

test:
	sbatch ${jobscript} test 2 1

timer:
	for cnm in 1000 2000 4000 6000 8000 10000 12000; do \
      sbatch simple_timer.sh $${cnm} ;\
  done

others_timer:
	for mym in ${others_of_timer}; do \
		for cnm in 1000 2000 4000 6000 8000 10000 12000; do \
				sbatch ${other_timer_script} $${mym} $${cnm} ;\
		done; \
   done

single_other_timer:
	sbatch ${other_timer_script} ${method} ${cnm}

single_other_series_timer:
	for cnm in 1000 2000 4000 6000 8000 10000 12000; do \
		 sbatch ${other_timer_script} ${method} $${cnm} ;\
  done


test_others_timer:
	for mym in ${others_of_timer}; do \
		for cnm in 1000; do \
				sbatch ${other_timer_script} $${mym} $${cnm} ;\
		done; \
   done

clean:
	rm job_simple*.err
	rm job_simple*.out
