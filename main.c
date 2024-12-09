/*  This is the main file of the interval selection experiments. This was used to
    generate the experiments described in the paper:

    "Interval Selection with Binary Predictions"
    
    Copyright (C) 2024  Christodoulos Karavasilis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.*/

#include <stdio.h>
#include <stdlib.h>
#include "algs.h"
#include <time.h>

#define ERROR_STEP 1000 // how many prediction elements to change per loop



int main(){
    int i, j, **t_arr, **inst,*seq, *opt, **perms, *sorted, opt_val;
    int *err_vals, *predictions, err_flag, step, ALG;
    int naive_max, naive_min, grNR_max, grNR_min, tk_max, tk_min, uPR_max, uPR_min, garay_min, garay_max, pPR_min, pPR_max, gar_mod_min, gar_mod_max, pPR_TWO_min, pPR_TWO_max;
    double naive_avg, grNR_avg, tk_avg, threeK_avg, uPR_avg, garay_avg, pPR_avg, gar_mod_avg, pPR_TWO_avg;
    FILE *naive_out, *grNR_out, *tk_out, *threeK_out, *uPR_out, *garay_out, *pPR_out, *garay_mod_out, *pPR_TWO_out;
    int cur_err, threeK_min, threeK_max;

    long curtime;
    curtime = time(NULL);
    srand((unsigned int) curtime);

    inst = malloc(N * sizeof(int *));
    for(i=0; i<N; i++) inst[i] = malloc(2 * sizeof(int));
    for(i=0; i<N; i++) inst[i][0] = -1;
    err_vals = malloc(N*sizeof(int)); arr_init(err_vals);
    predictions = malloc(N*sizeof(int)); arr_init(predictions);
    perms = malloc(PERM_NUM * sizeof(int *));
    for(i=0; i<PERM_NUM; i++) perms[i] = malloc(N*sizeof(int));


    seq = malloc(N*sizeof(int));
    opt = malloc(N*sizeof(int));
    arr_init(seq);arr_init(opt);


//////////////////////////////////////////////

    FILE *file = fopen("datasets/nasa_ipsc.txt", "r");
    if (file == NULL) {
        printf("Error opening file");
        return 1;
    }

    char line[MAX_LINE_LENGTH];
    int line_count = 0;

    while (fgets(line, sizeof(line), file) != NULL && line_count < MAX_LINES) {
        int num1, num2;
        if (sscanf(line, "%*d\t%d\t%*d\t%d", &num1, &num2) == 2) {
            if(num2 <= 0) continue;
            inst[line_count][0] = num1;
            inst[line_count][1] = num1 + num2;
            line_count++;
        }
    }
    fclose(file);
    printf("line count: %d\n\n",line_count);
    /* compute OPT value (proportional weights) */
    opt_val = prop_opt(inst,opt,line_count); printf("OPT is %d\n",opt_val);


    /////////
    for(i=0; i<line_count; i++) seq[i] = i;
////////////////////////////////////////////
    /** ---- UNCOMMENT THIS FOR UNIT ---- **/
    /*unit_opt(inst,opt); //compute optimal solution
    comp_error_unit(inst,opt,err_vals); //compute error per element
    for(i = 0; i < N; i++) {//init predictions with opt (perfect predictions)
        if(inst[i][0] == (-1)) break;
        if(opt[i] == 1) predictions[i] = 1;
        else predictions[i] = 0;
    }*/

    comp_error_prop(inst,opt,err_vals); //compute error per element
    for(i = 0; i < N; i++) {//init predictions with opt (perfect predictions)
        if(inst[i][0] == (-1)) break;
        if(opt[i] == 1) predictions[i] = 1;
        else predictions[i] = 0;
    }


    generate_perms(perms,line_count);
    //naive_out = fopen("naive_out.txt", "w");
    garay_out = fopen("garay_out.txt", "w");
    //tk_out = fopen("tk_out.txt", "w");
    //uPR_out = fopen("uPR_out.txt", "w");
    pPR_out = fopen("pPR_out.txt", "w");
    //grNR_out = fopen("grNR_out.txt", "w");
    //threeK_out = fopen("threeK_out.txt", "w");
    //pPR_TWO_out = fopen("pPR_two_out.txt", "w");

    /**************************************/
    /** Prediction-less algorithms first **/
    /**************************************/
    grNR_max = 0;
    grNR_min = N;
    grNR_avg = 0.0;
    tk_max = 0;
    tk_min = N;
    tk_avg = 0.0;
    garay_max = 0;
    garay_min = N;
    garay_avg = 0.0;
    gar_mod_avg = 0.0;
    gar_mod_min = N;
    gar_mod_max = 0;

    for(j = 0; j < PERM_NUM; j++){

            //ALG = greedy_noRev(inst,perms[j],0);
            //grNR_avg += (float)ALG / 10.0;
            //if (ALG > grNR_max) grNR_max = ALG;
            //if (ALG < grNR_min) grNR_min = ALG;

            ALG = garay_mod(inst,perms[j]);
            garay_avg += (double)ALG / 10.0;
            if (ALG > garay_max) garay_max = ALG;
            if (ALG < garay_min) garay_min = ALG;

            //ALG = garay_mod(inst,perms[j]);
            //gar_mod_avg += (float)ALG / 10.0;
            //if (ALG > gar_mod_max) gar_mod_max = ALG;
            //if (ALG < gar_mod_min) gar_mod_min = ALG;

            //ALG = two_k_alg(inst,perms[j]);
            //tk_avg += (float)ALG / 10.0;
            //if (ALG > tk_max) tk_max = ALG;
            //if (ALG < tk_min) tk_min = ALG;

    }
    //fprintf(tk_out, " %d,%d,%d,%f\n",cur_err,tk_max,tk_min,tk_avg);
    fprintf(garay_out, " %d,%d,%d,%f\n",cur_err,garay_max,garay_min,garay_avg);
    //fprintf(garay_mod_out, " %d,%d,%d,%f\n",cur_err,gar_mod_max,gar_mod_min,gar_mod_avg);
    //fprintf(grNR_out, " %d,%d,%d,%f\n",cur_err,grNR_max,grNR_min,grNR_avg);
    /*************************************/
    err_flag = 0;
    cur_err = 0;
    step = 0;
    /** Main Loop (over error values) **/
    while(err_flag == 0){
        printf("Step: %d\n",step);
        //////////////////////////////////
        /** init  values **/
        uPR_min = N;
        uPR_max = 0;
        uPR_avg = 0.0;
        threeK_min = N;
        threeK_max = 0;
        threeK_avg = 0.0;

        pPR_min = N;
        pPR_max = 0;
        pPR_avg = 0.0;
        pPR_TWO_min = N;
        pPR_TWO_max = 0;
        pPR_TWO_avg = 0.0;

        naive_min = N;
        naive_max = 0;
        naive_avg = 0.0;

        for(j = 0; j < PERM_NUM; j++){

            //ALG = naive_pr_algo(inst,perms[j],predictions,0);
            //naive_avg += (float)ALG / 10.0;
            //if (ALG > naive_max) naive_max = ALG;
            //if (ALG < naive_min) naive_min = ALG;

            ALG = prop_pr_rev(inst,perms[j],predictions);
            pPR_avg += (double)ALG / 10.0;
            if (ALG > pPR_max) pPR_max = ALG;
            if (ALG < pPR_min) pPR_min = ALG;

            //ALG = prop_pr_rev_TWO(inst,perms[j],predictions);
            //pPR_TWO_avg += (double)ALG / 10.0;
            //if (ALG > pPR_TWO_max) pPR_TWO_max = ALG;
            //if (ALG < pPR_TWO_min) pPR_TWO_min = ALG;


            //ALG = unit_pr_rev(inst,perms[j],predictions, opt, 1);
            //uPR_avg += (float)ALG / 10.0;
            //if (ALG > uPR_max) uPR_max = ALG;
            //if (ALG < uPR_min) uPR_min = ALG;

            //ALG = unit_pr_rev(inst,perms[j],predictions, opt, 0);
            //threeK_avg += (float)ALG / 10.0;
            //if (ALG > threeK_max) threeK_max = ALG;
            //if (ALG < threeK_min) threeK_min = ALG;

        }
        //fprintf(naive_out, " %d,%d,%d,%f\n",cur_err,naive_max,naive_min,naive_avg);
        //fprintf(uPR_out, " %d,%d,%d,%f\n",cur_err,uPR_max,uPR_min,uPR_avg);
        //fprintf(threeK_out, " %d,%d,%d,%f\n",cur_err,threeK_max,threeK_min,threeK_avg);
        fprintf(pPR_out, " %d,%d,%d,%f\n",cur_err,pPR_max,pPR_min,pPR_avg);
        //fprintf(pPR_TWO_out, " %d,%d,%d,%f\n",cur_err,pPR_TWO_max,pPR_TWO_min,pPR_TWO_avg);


        /////////////////////////////////////
        if( (step * ERROR_STEP) >=  line_count) err_flag = 1;//stop condition

        /** update error and predictions **/
        for(j = (step*ERROR_STEP); j < (step+1)*ERROR_STEP; j++){
            if( (predictions[j] == (-1)) || (j >= N) ){ //safeguard...
               break;
            }
            predictions[j] = 1 - predictions[j]; //swap prediction
            cur_err += err_vals[j];
        }
        step++;
    }

    //fclose(naive_out);
    fclose(garay_out);
    //fclose(tk_out);
    //fclose(threeK_out);
    //fclose(uPR_out);
    fclose(pPR_out);
    //fclose(pPR_TWO_out);
    //fclose(grNR_out);
    //fclose(garay_mod_out);

    for(i=0; i<N; i++) free(inst[i]);
    free(inst);
    free(seq); free(opt); free(err_vals); free(predictions);
    for(i=0; i<PERM_NUM; i++) free(perms[i]);
    free(perms);

////////////////////////
    return 0;
}
