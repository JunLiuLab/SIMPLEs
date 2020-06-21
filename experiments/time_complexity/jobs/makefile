other_timer_script := others_timer.sh
others_of_timer := rmagic saver scimpute scrabble viper
method ?= control
cnm ?= 12000

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
