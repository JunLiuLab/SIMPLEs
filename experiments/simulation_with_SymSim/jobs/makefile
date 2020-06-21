rscript := eval_simulation_singlejob.R
jobscript := eval_simulation_singlejob.sh
other_methods := control scrabble viper saver rmagic scimpute
method ?= control
myrep ?= 20

test:
	-mkdir result
	Rscript --vanilla ../${rscript} -m ${method} --nrep 2

rmagic:
	-mkdir result
	Rscript --vanilla ../${rscript} -m rmagic --nrep ${myrep}

tsymsim:
	Rscript --vanilla -e "library(SymSim)"

tsub:
	sbatch ${jobscript} ${method} 2 1 2

sub_other_methods:
	$(foreach var,$(other_methods), sbatch ${jobscript} $(var) 2 1 ${myrep};)

sub_simples:
	for k in 2 4 8 10; do \
     for m in 1 2 4 5; do \
        sbatch ${jobscript} simple $${k} $${m} ${myrep};\
     done; \
  done

tsub_simples:
	for k in 2 4; do \
     for m in 1; do \
        sbatch ${jobscript} simple $${k} $${m} ${myrep};\
     done; \
  done

scvi:
	Rscript --vanilla ../${rscript} -m scvi --nrep ${myrep}

del_rs:
	rm -rf result/*

clean:
	rm job_simple*.err
	rm job_simple*.out
