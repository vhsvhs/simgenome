#include "common.h"

int main( int argc, char **argv )
{
	srand(time(0));

	print_splash();

	/* Read the command-line and load user-specified files */
	settings *ss = make_settings();
	read_cli(argc, argv, ss);
	print_settings(ss);

	/* There are three ways to make a population:
	 * 1. Load a saved population using the deserialization method.
	 * 2. Initialize a random population.
	 * 3. Default: Build a single genome by reading genes from a file,
	 * then clone that individual N times to make the population.
	 */
	t_pop *pop;
	if (ss->load_save_pop == true){
		pop = deserialize_population(ss);
	}
	else if(ss->build_random_population == true){
		if (ss->verbosity > 0){	 printf("\n. I'm building a random population....\n"); }
		t_genome *gn = make_genome_random(ss);
		pop = make_population(gn, ss);
	}
	else{
		int ngenes;
		t_gene** mygenes = read_genes_from_file(ss, ngenes);
		if (ss->verbosity > 0){	printf("\n. I found %d genes.\n", ngenes); }

		t_genome *gn = make_genome(ngenes, mygenes, ss);
		if (ss->verbosity > 0){	 printf("\n. I'm building a population....\n"); }

		pop = make_population(gn, ss);
	}


	if (ss->verbosity > 0){	 printf("\n\n. I'm building a landscape....\n"); }
	t_landscape *l = (t_landscape *)malloc(sizeof(t_landscape));

	t_ruleset** rs = read_rulesets_from_file(ss, l->nrulesets, l->ntime);
	l->rulesets = rs;
	if (ss->verbosity > 0){	 printf("\n. I found %d rulesets.\n", l->nrulesets); }

	if (ss->verbosity > 0){	 printf("\n. I'm building the simulation. . .\n");  }
	t_ga* ga = make_ga();
	ga->pop = pop;
	ga->l = l;

	if ( !is_world_consistent(ga) ){
		printf("\n. Something is inconsistent with your setup.\n");
	}

	if (ss->verbosity > 0){	 printf("\n. The simulation is starting now.\n"); }
	runsim(ga, ss);

	if (ss->verbosity > 0){	 printf("\n. The simulation is finished."); }

	if (ss->verbosity > 0){	 printf("\n. Goodbye\n\n"); }

	free_ga( ga );
	free_settings(ss);

	return 1;
}
