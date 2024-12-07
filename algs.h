
#define MAX_LINE_LENGTH 200
#define MAX_LINES 100000
#define PERM_NUM 10
#define N 77500



/**** No Revoking ****/
//naive_pr_algo             //Same for both unit and prop //w_type = 0 for unit, 1 for prop
int naive_pr_algo(int **inst, int *seq, int *pred, int w_type);

//greedy_noRev             //Same for both unit and prop
int greedy_noRev(int **inst, int *seq, int w_type);

/**** Revoking ****/
//unit_pr_rev              // (predictions) PAPER
int unit_pr_rev(int **inst, int *seq, int *pred, int *opt, int m_flag);//m_flag = 1 means marks carries over

//prop_pr_rev              // (predictions) PAPER
int prop_pr_rev(int **inst, int *seq, int *pred);
int prop_pr_rev_TWO(int **inst, int *seq, int *pred);//different choice of lambda

//unit_rev                  // 2k algo
int two_k_alg(int **inst, int *seq);

//prop_rev                  // Garay
int garay(int **inst, int *seq);
int garay_mod(int **inst, int *seq); //modifed Garay comparing to total conflict


/**** OPTs ****/
//unit_opt              //Faigle
void unit_opt(int **inst, int *opt);

//prop_opt              //Tardos DP
int prop_opt(int **inst, int *opt, int lines); //also returns objective value of opt solution


/**** HELPER FUNCTIONS ****/

int conf_left(int **inst, int cur, int *sol, int *c); //// 0 ==> no conflict, 1 ==> partial conflict, 2==> containment
int conf_right(int **inst, int cur, int *sol, int *c); //// 0 ==> no conflict, 1 ==> partial conflict, 2==> containment
int conf_any(int **inst, int cur, int *sol); //0 ==> no conflict
int conf_rightmost(int **inst, int cur, int *sol, int *c);//0 ==> no conflict, else 1

void find_conflicts(int **inst, int cur, int *sol, int *leftmost, int *rightmost); //store leftmost and rightmost conflicting intervals

//compute objective value given a solution (need one for unit and one for prop)
int obj_val_unit(int *sol);
int obj_val_prop(int **inst, int *sol);

//init array with -1
void arr_init(int *ar);

//compute the error of each element, given a fixed (optimal) solution
void comp_error_unit(int **inst, int *sol, int *err_vals); /* sol should be an optimal solution (which can be interpreted as predictions) */
void comp_error_prop(int **inst, int *sol, int *err_vals);

void generate_perms(int **perms, int n);//n <= N is actual number of items. perms is PERM_NUMxN

////////////////////////////
int compare(const void *a, const void *b); //used to quicksort array
void get_sorted_arr(int **inst, int lines, int *sorted);
//the following follows notation from the DP algorithm ('Algorithm Design' Kleinberg & Tardos)
void compute_pjs(int **inst, int *sorted, int lines, int *pj_array);
void compute_M_array(int **inst, int *sorted, int lines, int *pj_array, int *M_array);
