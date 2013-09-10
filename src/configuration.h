#define DEF_VERBOSITY 2

#define N_STATES 4 /* for nucleotide regulatory sequences */

#define GENE_NAME_MAX 100
#define FILEPATH_LEN_MAX 300

/* MAXLEN is the max line length for user-specified file.
 * This value should always be longer than URSLEN_DEFAULT
 */
#define MAXLEN 1200
#define MAX_TOKENS 6
#define MAX_PSAM_LEN 100

#define MAX_DDG 2.0
#define MIN_DDG -2.0

#define PWMLENMU 0.3
#define PWMLENMUMAX 2
#define DDGMU 0.5

#define URSLEN_DEFAULT 100
#define PSAMLEN_DEFAULT 6
#define NGENES_DEFAULT 10

/* Output directory path defaults */
#define POPDIR "LOG_POPS"
#define EXPRDIR "LOG_EXPR"
#define FITNESS_DIR "LOG_FITNESS"
#define COOP_DIR "LOG_COOP"
#define DBD_DIR "LOG_DBD"
#define CONFIG_DIR "LOG_OCCUPANCY"
