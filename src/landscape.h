/* A t_rule struct defines a single fitness rule that must
 * be met in order for maximum fitness.
 */
typedef struct __Rule {
	int timepoint;
	int repid; /* The reporter ID */
	double expr_level;
	int rule_type; /* 0 = greater than or equal to, 1 = less than or equal to */
	double weight; /* Weight this rule differently from other rules when computing fitness. Default = 1.0 */
}t_rule;

typedef struct __Input{
	int start; 		/* start timepoint */
	int stop; 		/* stop timepoint */
	int gid; 		/* The gene ID */
	double expr_level; /* expression level of gid */
}t_input;

/* t_ruleset is a collection of rules that must be simultaneously
 * satisfied in order to achieve maximum fitness.
 */
typedef struct __Ruleset{
	int id; /* ID of this ruleset */
	t_input** inputs;
	int ninputs;
	t_rule** rules;
	int nrules;
}t_ruleset;

typedef struct __Landscape {
	t_ruleset** rulesets;
	int nrulesets;
	int ntime;
	//int* r; /* r[gene id] = PSAM length */
}t_landscape;


t_rule* make_rule(int timepoint, int repid, double expr, int ruletype, double weight);
t_input* make_input(int start, int stop, int gid, double expr);
t_ruleset* make_ruleset(int id, int nrules, int ninputs);
t_ruleset** read_rulesets_from_file(settings* ss, int &ret_n, int &ntime);
void free_landscape(t_landscape* l);
