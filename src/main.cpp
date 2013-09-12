#include "common.h"

int main( int argc, char **argv )
{
	srand(time(0));

	settings *ss = make_settings();
	read_cli(argc, argv, ss);
	print_settings(ss);
	exit(1);

	char urspath[] = "/Users/victor/Applications/simgenome-c-beta/examples/urs.time092013.txt";
	ss->urspath = urspath;

	char psampath[] = "/Users/victor/Applications/simgenome-c-beta/examples/five.psam";
	ss->psampath = psampath;

	int ngenes;
	t_gene** mygenes = read_genes_from_file(ss, ngenes);
	printf("\n. I found %d genes.\n", ngenes);

	t_genome *gn = make_genome(ngenes, mygenes, ss);
	printf("\n. I'm building a population....\n");
	t_pop *pop = make_population(1, gn, ss);

	printf("\n. I'm building a landscape....\n");
	t_landscape *l = (t_landscape *)malloc(sizeof(t_landscape));

	char rulepath[] = "/Users/victor/Applications/simgenome-c-beta/examples/five.rules.txt";
	ss->rulepath = rulepath;
	t_ruleset** rs = read_rulesets_from_file(ss, l->nrulesets, l->ntime);
	l->rulesets = rs;
	printf("\n. I found %d rulesets.\n", l->nrulesets);

	printf("\n. I'm building a G.A.\n");
	t_ga* ga = make_ga();
	ga->pop = pop;
	ga->l = l;

	printf("\n. I'm running the G.A.\n");
	runsim(ga, ss);

	printf("\n. I'm freeing the population....\n");
	free_pop( pop );

	return 1;
}
