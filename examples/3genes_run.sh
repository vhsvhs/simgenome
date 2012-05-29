mpirun -np 2 python runme.py --patternpath ./examples/3genes_rules.txt --verbose 8 --runid 3genes.test --growth_rate 6.0 --decay_rate 4.0 --popsize 6 --mu 0.1 --elitemu 0.0 --pwmmu 0.0 --maxtime 10 --eliteproportion 0.2 --numtr 2 --numreporter 1 --init_pwm_len 4 --init_urs_len 10 --iid_samples 10 --pop_path 3genes.test/POP_HISTORY/population.gen4.pickle --start_generation 5