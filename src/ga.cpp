#include "common.h"

t_ga* make_ga(){
	t_ga* ga;
	ga = (t_ga*)malloc(sizeof(t_ga));
	return ga;
}

/* Assumes that ga has already been built and
 * that ga->l and ga->pop exist.
 */
void runsim(t_ga* ga, settings* ss){

	/* For each generation */
	int start_gen = ss->gen_counter;
	for (int ii = start_gen;
			ii < (start_gen + ss->max_gens);
			ii++){

		ss->gen_counter = ii;

		/*
		 * Pre-fitness Setup
		 */
		if (ii == 0){ // Only do this for generation zero
			for(int ii = 0; ii < ga->pop->ngenomes; ii++){
				// allocate memory AND init gene expression
				build_lifespan(ga->pop->genomes[ii], ga->l->ntime);
			}
		}
		else{
			for (int gid = 0; gid < ga->pop->ngenomes; gid++){
				// just re-init gene expression (don't re-allocate memory)
				reset_lifespan( ga->pop->genomes[gid]);
			}
		}

		/*
		 * FITNESS
		 */

		/* f is an array of fitness values, one for each individual */
		double* f = (double *)malloc(ga->pop->ngenomes * sizeof(double));
		for (int qq = 0; qq < ga->pop->ngenomes; qq++){
			f[qq] = 0.0;
		}
		double maxf = 0;
		double minf = 1.0;
		double meanf = 0;
		double medianf = 0;
		double stdf = 0;

		/* Consider each individual */
		for (int gid=0; gid < ga->pop->ngenomes; gid++){
			/*
			 * GET FITNESS of an individual
			 */
			f[gid] = get_fitness(ga->pop->genomes[gid], ga->l, ss);

			// to-do: write the expression cran
		}

		printf("ga 65\n");

		get_fitness_stats( f, ga->pop->ngenomes,
				maxf, minf, meanf, medianf, stdf);

		printf("==================================\n");
		printf(". Generation %d Fitness:\n", ii);
		for (int gid = 0; gid < ga->pop->ngenomes; gid++){
			printf("\t ID %d, f= %f\n", gid, f[gid]);
		}
		printf("==================================\n");

		log_fitness(f, ga->pop->ngenomes, ss);

		exit(1);

		/*
		 * SELECTIVELY REPRODUCE
		 */
		reproduce( ga->pop, ss, f);

		/*
		 * MUTATION
		 */
		mutate(ga->pop, ss);
	}
}

void get_fitness_stats(double *f, int len, double &max, double &min, double &mean, double &median, double &stddev) {
	max = f[0];
	min = f[0];
	double sum = 0.0;
	for (int ii = 0; ii < len; ii++){
		if (f[ii] > max) { max = f[ii]; }
		else if (f[ii] < min) { min = f[ii]; }
		sum += f[ii];
	}
	mean = sum / len;
	// to-do: median and stdev
}

void free_ga(t_ga* ga){
	free_pop(ga->pop);
	free_landscape(ga->l);
}
