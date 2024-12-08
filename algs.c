#include "algs.h"
#include <stdlib.h>


/** No Revoking **/
//naive_pr_algo             //Same for both unit and prop
int naive_pr_algo(int **inst, int *seq, int *pred, int w_type){
    int i, *sol, ALG, cur_item;

    sol = malloc(N*sizeof(int));
    arr_init(sol);
    ////////////////////
    for(i=0; i<N; i++){
        cur_item = seq[i];
        if(cur_item == (-1)) break;
        if(conf_any(inst,cur_item,sol) == 1) continue; //if any conflict --> skip
        if(pred[cur_item] == 1) sol[cur_item] = 1; //take it if predicted to be optimal
    }
    //////////////////////////
    if(w_type == 0)
        ALG = obj_val_unit(sol);
    else
        ALG = obj_val_prop(inst,sol);
    free(sol);
    return ALG;
    //////////////
}

//greedy_noRev             //Same for both unit and prop
int greedy_noRev(int **inst, int *seq, int w_type){
    int i, *sol, ALG, cur_item;

    sol = malloc(N*sizeof(int));
    arr_init(sol);
    //////////////////////
    for(i=0; i<N; i++){
        cur_item = seq[i];
        if(cur_item == (-1)) break;
        if(conf_any(inst,cur_item,sol) == 1) continue; //if any conflict --> skip
        sol[cur_item] = 1; //take it if no conflict
    }
    //////////////////////////
    if(w_type == 0)
        ALG = obj_val_unit(sol);
    else
        ALG = obj_val_prop(inst,sol);
    free(sol);
    return ALG;
    //////////////
}



/** Revoking **/

//unit_pr_rev              // (predictions) PAPER
int unit_pr_rev(int **inst, int *seq, int *pred, int *opt, int m_flag){
    int i, *sol, *marked, ALG, cur_item, conf_item_1, conf_item_2, l_c, r_c;

    sol = malloc(N*sizeof(int));
    marked = malloc(N*sizeof(int));
    arr_init(sol);arr_init(marked);
    /////////////////////////////////
    for(i=0; i<N; i++){
        cur_item = seq[i];
        if(cur_item == (-1)) break;
        l_c = conf_left(inst,cur_item,sol,&conf_item_1);
        r_c = conf_right(inst,cur_item,sol,&conf_item_2);
        if(l_c == 2){ //containment --> replace
            sol[conf_item_1] = -1;
            sol[cur_item] = 1;
            if( (marked[conf_item_1] == 1) && (m_flag == 1) ) marked[cur_item] = 1;//carry the mark to new interval
        } else if((l_c == 0) && (r_c == 0)){ //no conflict --> simply accept
            sol[cur_item] = 1;
        } else if(r_c != 2){ //only partial conflicts
            if(pred[cur_item] != 1) continue;//if not predicted to be optimal, skip
            if(l_c != 0){//there's a partial left conflict
                if(marked[conf_item_1] == 1) //continue if marked
                    continue;
            }
            if(r_c != 0){//there's a partial right conflict
                if(marked[conf_item_2] == 1) //continue if marked
                    continue;
            }
            //if none of the partial conflicts are marked, we can replace them
            if (l_c != 0) sol[conf_item_1] = -1;
            if (r_c != 0) sol[conf_item_2] = -1;
            //accept new item and mark it
            sol[cur_item] = 1;
            marked[cur_item] = 1;
        }
    }
    /////////////////////////
    ALG = obj_val_unit(sol);
    free(sol);free(marked);
    return ALG;
    //////////////
}

//prop_pr_rev              // (predictions) PAPER
int prop_pr_rev(int **inst, int *seq, int *pred){
    int i, j, *sol, ALG, cur_item, conf_item_L, conf_item_R, w_c, w_I;
    double lambda, rho;
    lambda = 1.618; //choice of lambda
    rho = 0.5; //make it rho-increasing

    sol = malloc(N*sizeof(int));
    arr_init(sol);
    /////////////////////////////////
    for(i=0; i<N; i++){
        cur_item = seq[i];
        if(cur_item == (-1)) break;

        w_c = 0; //init total conflicting weight
        w_I = inst[cur_item][1] - inst[cur_item][0]; //weight of new interval
        find_conflicts(inst,cur_item,sol,&conf_item_L,&conf_item_R);
        for(j = conf_item_L; j <= conf_item_R; j++){//loop over all conflicting items
            if(sol[j] != 1) continue;
            w_c += inst[j][1] - inst[j][0];
        }
        if( (double)w_I >= lambda * (double)w_c){//main replace rule
            for(j = conf_item_L; j <= conf_item_R; j++) sol[j] = -1;//replace conflicts
            sol[cur_item] = 1; //accept new item
        } else if(pred[cur_item] == 1){//fails main rule but considered optimal by predictions
            if( (double)w_I < rho * (double)w_c ) continue; //fails even the relaxed rule
            for(j = conf_item_L; j <= conf_item_R; j++){
                if ( (pred[j] == 1) && (sol[j] == 1) ) continue; //conflicts with predicted optimal interval in solution
            }
            for(j = conf_item_L; j <= conf_item_R; j++) sol[j] = -1;//replace conflicts
            sol[cur_item] = 1; //accept new item
        }
    }
    /////////////////////////////
    ALG = obj_val_prop(inst,sol);
    free(sol);
    return ALG;
    //////////////
}

int prop_pr_rev_TWO(int **inst, int *seq, int *pred){
    int i, j, *sol, ALG, cur_item, conf_item_L, conf_item_R, w_c, w_I;
    double lambda;
    lambda = 4.0; //choice of lambda

    sol = malloc(N*sizeof(int));
    arr_init(sol);
    /////////////////////////////////
    for(i=0; i<N; i++){
        cur_item = seq[i];
        if(cur_item == (-1)) break;

        w_c = 0; //init total conflicting weight
        w_I = inst[cur_item][1] - inst[cur_item][0]; //weight of new interval
        find_conflicts(inst,cur_item,sol,&conf_item_L,&conf_item_R);
        for(j = conf_item_L; j <= conf_item_R; j++){//loop over all conflicting items
            if(sol[j] != 1) continue;
            w_c += inst[j][1] - inst[j][0];
        }
        if( (double)w_I >= lambda * (double)w_c){//main replace rule
            for(j = conf_item_L; j <= conf_item_R; j++) sol[j] = -1;//replace conflicts
            sol[cur_item] = 1; //accept new item
        } else if(pred[cur_item] == 1){//fails main rule but considered optimal by predictions
            if( w_I < w_c ) continue; //fails even the relaxed rule
            for(j = conf_item_L; j <= conf_item_R; j++){
                if ( (pred[j] == 1) && (sol[j] == 1) ) continue; //conflicts with predicted optimal interval in solution
            }
            for(j = conf_item_L; j <= conf_item_R; j++) sol[j] = -1;//replace conflicts
            sol[cur_item] = 1; //accept new item
        }
    }
    /////////////////////////////
    ALG = obj_val_prop(inst,sol);
    free(sol);
    return ALG;
    //////////////
}

int garay_mod(int **inst, int *seq){
        int i, j, *sol, ALG, cur_item, conf_item_L, conf_item_R, w_c, w_I;
    double lambda;
    lambda = 1.0; //choice of lambda

    sol = malloc(N*sizeof(int));
    arr_init(sol);
    /////////////////////////////////
    for(i=0; i<N; i++){
        cur_item = seq[i];
        if(cur_item == (-1)) break;

        w_c = 0; //init total conflicting weight
        w_I = inst[cur_item][1] - inst[cur_item][0]; //weight of new interval
        find_conflicts(inst,cur_item,sol,&conf_item_L,&conf_item_R);
        for(j = conf_item_L; j <= conf_item_R; j++){//loop over all conflicting items
            if(sol[j] != 1) continue;
            w_c += inst[j][1] - inst[j][0];
        }
        if( (double)w_I >= lambda * (double)w_c){//main replace rule
            for(j = conf_item_L; j <= conf_item_R; j++) sol[j] = -1;//replace conflicts
            sol[cur_item] = 1; //accept new item
        }
    }
    /////////////////////////////
    ALG = obj_val_prop(inst,sol);
    free(sol);
    return ALG;
    //////////////
}

//unit_rev                  // 2k algo
int two_k_alg(int **inst, int *seq){
    int i, *sol, ALG, cur_item, conflicting_item, l_c;

    sol = malloc(N*sizeof(int));
    arr_init(sol);
    /////////////////////////////////
    for(i=0; i<N; i++){
        cur_item = seq[i];
        if(cur_item == (-1)) break;

        if(conf_right(inst,cur_item,sol,&conflicting_item) != 0) continue;
        l_c = conf_left(inst,cur_item,sol,&conflicting_item);
        if(l_c == 1) continue; //if partial conflict, skip
        else if(l_c == 2){ //containment --> replace
            sol[conflicting_item] = -1; //remove conflicting interval
            sol[cur_item] = 1;
        } else sol[cur_item] = 1;//no conflict --> simply accept
    }
    /////////////////////////////
    ALG = obj_val_unit(sol);
    free(sol);
    return ALG;
    //////////////
}

//prop_rev                  // Garay
int garay(int **inst, int *seq){
    int i, j, *sol, ALG, cur_item, conf_item_L, conf_item_R, w_max;

    sol = malloc(N*sizeof(int));
    arr_init(sol);
    /////////////////////////////////
    for(i=0; i<N; i++){
        cur_item = seq[i];
        if(cur_item == (-1)) break;

        w_max = 0; //init max conflicting weight
        find_conflicts(inst,cur_item,sol,&conf_item_L,&conf_item_R);
        for(j = conf_item_L; j <= conf_item_R; j++){//loop over all conflicting items
            if(sol[j] != 1) continue;
            if( (inst[j][1] - inst[j][0]) > w_max) w_max = inst[j][1] - inst[j][0];
        }
        if( (double)(inst[cur_item][1] - inst[cur_item][0]) > 1.618 * (double)w_max){//Garay replace rule
            for(j = conf_item_L; j <= conf_item_R; j++) sol[j] = 0;//replace conflicts
            sol[cur_item] = 1; //accept new item
        }
    }
    /////////////////////////////
    ALG = obj_val_prop(inst,sol);
    free(sol);
    return ALG;
    //////////////
}




/** OPTs **/
//unit_opt              //Faigle
void unit_opt(int **inst, int *opt){
    int i, cur_job;

    cur_job = 0;
    opt[0] = 1; //start with first interval in the solution
    for(i = 1; i < N; i++){
        if (inst[i][0] == -1) break; //stop if instance is over

        if(inst[i][0] >= inst[cur_job][1]){//new interval begins after current one ends
            cur_job = i;//add new interval to the solution
            opt[i] = 1;
        } else if (inst[i][1] < inst[cur_job][1]){//contained in current job ==> replace
            opt[cur_job] = -1;
            cur_job = i;
            opt[cur_job] = 1;
        }
    }
}

//prop_opt              //Tardos DP
int prop_opt(int **inst, int *opt, int lines){
    int i, *sorted, *pj_array, *M_array, *tmp_opt, opt_val;

    sorted = malloc(N*(sizeof(int))); arr_init(sorted);
    pj_array = malloc(N*(sizeof(int))); arr_init(pj_array);
    M_array = malloc(N*(sizeof(int))); arr_init(M_array);
    tmp_opt = malloc(N*(sizeof(int))); arr_init(tmp_opt);

    get_sorted_arr(inst,lines,sorted);//computed sorted array
    compute_pjs(inst,sorted,lines,pj_array);
    compute_M_array(inst,sorted,lines,pj_array,M_array);
    ///////////////
    /** Compute OPT solution given M **/
    i = lines-1;
    while(i >= 0){
        if(i == 0){
            tmp_opt[0] = 1;
            i--;
            break;
        }
        int m_p_j;
        if(pj_array[i] >= 0) m_p_j = M_array[pj_array[i]];
        else m_p_j = 0;
        
        if( (inst[ sorted[i] ][1] - inst[ sorted[i] ][0]) + m_p_j >= M_array[i-1]){
            //this element is optimal, move to its pj
            tmp_opt[i] = 1;
            i = pj_array[i];
        } else i--;//not optimal element, move to previous onee
    }
    ////////////////
    opt_val = M_array[lines-1]; 
    for(i=0; i<lines; i++){ //updated main solution
        if(tmp_opt[i] == 1) opt[sorted[i]] = 1;
    }
    free(sorted);free(pj_array);free(M_array);free(tmp_opt);
    /////////////////
    return opt_val;
}




/** HELPER FUNCTIONS **/

int conf_left(int **inst, int cur, int *sol, int *c){// 0 ==> no conflict, 1 ==> partial conflict, 2==> containment
    int i, conflict;
    conflict = 0;
    for(i = cur-1; i>=0; i--){
        if(sol[i] != 1) continue; //if the interval is not in the current solution, skip

        if((inst[i][1] > inst[cur][0]) &&(inst[i][1] < inst[cur][1]) ){//partial left conflict
            conflict = 1;
            (*c) = i;
            break;
        }

        if(inst[i][1] >= inst[cur][1]){ //proper containment <==> i contains cur
            conflict = 2;
            (*c) = i;
            break;
        }
    }
    return conflict;
}


int conf_right(int **inst, int cur, int *sol, int *c){ // 0 ==> no conflict, 1 ==> partial conflict, 2==> containment
    int i, conflict;
    conflict = 0;
    for(i = cur+1; i<N; i++){
        if(inst[i][0] == (-1)) break;//instance is finished...
        if(inst[i][0] >= inst[cur][1]) break; //there cannot be conflict if interval starts after current one ends
        if(sol[i] != 1) continue; //if the interval is not in the current solution, skip

        if(inst[i][1] <= inst[cur][1]){//cur contains i
            conflict = 2;
            (*c) = i;
            break;
        }
        if(inst[i][1] > inst[cur][1]){//partial conflict
            conflict = 1;
            (*c) = i;
            break;
        }
    }
    return conflict;
}


int conf_any(int **inst, int cur, int *sol){ //0 ==> no conflict
    int c;
    if(conf_left(inst,cur,sol,&c) != 0)
        return 1;
    if(conf_right(inst,cur,sol,&c) != 0)
        return 1;
    return 0;
}

int conf_rightmost(int **inst, int cur, int *sol, int *c){
    int i, conflict;
    conflict = 0;
    for(i = cur+1; i<N; i++){
        if(inst[i][0] == (-1)) break;//instance is finished...
        if(inst[i][0] >= inst[cur][1]) break; //there cannot be conflict if interval starts after current one ends
        if(sol[i] != 1) continue; //if the interval is not in the current solution, skip

        if(inst[i][1] <= inst[cur][1]){//cur contains i
            conflict = 1; //just 1 instead of 2 (conf_right case), as we don't care about its exact nature
            (*c) = i;
            //break;
        }
        if(inst[i][1] > inst[cur][1]){//partial conflict
            conflict = 1;
            (*c) = i;
            //break;
        }
    }
    return conflict;
}


void find_conflicts(int **inst, int cur, int *sol, int *leftmost, int *rightmost){ //store leftmost and rightmost conflicting intervals
    int i, l_c, r_c;

    l_c = conf_left(inst,cur,sol,leftmost);
    if (l_c == 2){ //cur is contained in an interval in the solution
        (*rightmost) = (*leftmost);
        return;
    } else if (l_c == 1){ //partially conflicting interval on the left
        r_c = conf_rightmost(inst,cur,sol,rightmost);
        if(r_c == 0){
            (*rightmost) = cur; //hack...
            return;
        } else return;

    } else{ //no conflict on the left
        (*leftmost) = cur; //hack...
        r_c = conf_rightmost(inst,cur,sol,rightmost);
        if(r_c == 0){//no conflict on right
            (*rightmost) = cur; //hack...
            return;
        } else return; //rightmost is already set to something
    }
}



//compute objective value given a solution (need one for unit and one for prop)
int obj_val_unit(int *sol){
    int i, ALG;
    ALG = 0;
    for(i=0; i<N; i++){
        if(sol[i] == 1) ALG++; //add interval to ALG
    }
    return ALG;
}

int obj_val_prop(int **inst, int *sol){
    int i, ALG;
    ALG = 0;
    for(i=0; i<N; i++){
        if(sol[i] == 1) ALG += (inst[i][1] - inst[i][0]); //add length of interval to ALG
    }
    return ALG;
}

void arr_init(int *ar){
    int i;
    for(i=0; i<N; i++) ar[i] = -1;
}


void comp_error_unit(int **inst, int *sol, int *err_vals){
    int i, j, l_c, r_c;

    for(i=0; i<N; i++){
        if(inst[i][0] == (-1)) break; //we reached end of input
        if(sol[i] == 1) err_vals[i] = 1; //wrongly predicted to be nonoptimal
        else{ //wrongly predicted to be optimal
            find_conflicts(inst,i,sol,&l_c,&r_c);
            err_vals[i] = 0;
            for(j=l_c; j<=r_c; j++){
                if(sol[j] == 1) err_vals[i]++;
            }
            err_vals[i]--;//subtract its own weight
        }
    }
}

void comp_error_prop(int **inst, int *sol, int *err_vals){
    int i, j, l_c, r_c;

    for(i=0; i<N; i++){
        if(inst[i][0] == (-1)) break; //we reached end of input
        if(sol[i] == 1) err_vals[i] = inst[i][1] - inst[i][0]; //wrongly predicted to be nonoptimal
        else{ //wrongly predicted to be optimal
            find_conflicts(inst,i,sol,&l_c,&r_c);
            err_vals[i] = 0;
            for(j=l_c; j<=r_c; j++){
                if(sol[j] == 1) err_vals[i] += inst[j][1] - inst[j][0];
            }
            err_vals[i] -= inst[i][1] - inst[i][0]; //subtract its own weight
        }
    }
}


void generate_perms(int **perms, int n){
    int i, j, item, temp;

    for(i = 0; i < PERM_NUM; i++){
        for(j = 0; j < n; j++){
            perms[i][j] = j;//initialize all permutations 0,1,2,...n
        }
    }
    /* iterate over everything again, picking a random value out of the remaining array */
    for(j = 0; j < n; j++){
        for(i = 0; i < PERM_NUM; i++){
            item = j + rand()%(n-j);
            temp = perms[i][j];
            perms[i][j] = perms[i][item];
            perms[i][item] = temp;
        }
    }
}

int compare(const void *a, const void *b) {
    const int *row_a = *(const int **)a;
    const int *row_b = *(const int **)b;

    return row_a[1] - row_b[1];
}


void get_sorted_arr(int **inst, int lines, int *sorted){
    int i, **temp_arr;

    temp_arr = malloc(N*sizeof(int *));
    for(i=0; i<N; i++) temp_arr[i] = malloc(2*sizeof(int));
    for(i=0; i<N; i++){
        temp_arr[i][0] = i;//inst[i][0];
        temp_arr[i][1] = inst[i][1];
    }

    /////////////////////////////////////
    //sort
    qsort(temp_arr, lines, sizeof(int) *2, compare);

    /////////////////////////////

    for(i=0; i<lines; i++) sorted[i] = temp_arr[i][0];

    /////////////////////////////////
    for(i=0; i<N; i++) free(temp_arr[i]);
    free(temp_arr);
    ///////////////////////////////////////////////
}



void compute_pjs(int **inst, int *sorted, int lines, int *pj_array){
    int i, temp, pj_val;

    for(i=lines-1; i>=0; i--){
        pj_val = -1;
        temp = i-1;
        while( (temp >= 0) && (inst[ sorted[temp] ][1] > inst[ sorted[i] ][0] ) )
            temp--;
        if(temp >= 0){ //there exists an interval that can be scheduled...
            pj_val = temp;
        }
        pj_array[i] = pj_val;//-1 if all previous ones conflict, else the index
    }
}

void compute_M_array(int **inst, int *sorted, int lines, int *pj_array, int *M_array){
    int i, val_a, val_b;

    for(i = 0; i < lines; i++){
        /* compute val_a */
        if(pj_array[i] != (-1)) val_a = (inst[sorted[i]][1] - inst[sorted[i]][0]) + M_array[pj_array[i]]; //if M[p(j)] > 0
        else val_a = inst[sorted[i]][1] - inst[sorted[i]][0];

        /* compute val_b */
        if(i > 0) val_b = M_array[i - 1];
        else val_b = 0;

        if(val_a > val_b) M_array[i] = val_a;
        else M_array[i] = val_b;
    }
    ///////////////////////
}




/////////////////////////////////////////////////////////










