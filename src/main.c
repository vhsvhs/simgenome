#include "common.h"

int main( int argc, char **argv )
{
	srand(time(0));

	print_splash();
	settings *ss = make_settings();
	read_cli(argc, argv, ss);
	print_settings(ss);

	int ngenes;
	t_gene** mygenes = read_genes_from_file(ss, ngenes);
	if (ss->verbosity > 0){	printf("\n. I found %d genes.\n", ngenes); }

	t_genome *gn = make_genome(ngenes, mygenes, ss);
	if (ss->verbosity > 0){	 printf("\n. I'm building a population....\n"); }


	t_pop *pop = make_population(gn, ss);

	if (ss->verbosity > 0){	 printf("\n. I'm building a landscape....\n"); }
	t_landscape *l = (t_landscape *)malloc(sizeof(t_landscape));

	t_ruleset** rs = read_rulesets_from_file(ss, l->nrulesets, l->ntime);
	l->rulesets = rs;
	if (ss->verbosity > 0){	 printf("\n. I found %d rulesets.\n", l->nrulesets); }

	if (ss->verbosity > 0){	 printf("\n. I'm building a G.A.\n");  }
	t_ga* ga = make_ga();
	ga->pop = pop;
	ga->l = l;

	if (ss->verbosity > 0){	 printf("\n. I'm running the G.A.\n"); }
	runsim(ga, ss);

	if (ss->verbosity > 0){	 printf("\n. The simulation is finished."); }
	free_ga( ga );
	if (ss->verbosity > 0){	 printf("\n. Goodbye\n\n"); }

	free_settings(ss);

	return 1;
}
