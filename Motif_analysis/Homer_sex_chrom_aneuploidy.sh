#!/bin/bash

# A script to run homer de novo motif analysis.

# Job submission function
submit_job() {
    local input_file="$1"
    local output_dir="$2"
    local bg_file="$3"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=homer
#SBATCH --output=homer_%j.out.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --partition=20    
#SBATCH --mail-type=END              

findMotifs.pl "$input_file" human "$output_dir" -bg "$bg_file" -p 40 -start -1000 -end 1000 -nogo -fdr 100
EOF
}

# Submit jobs
submit_job "lcl_x_resp_up_genes_auto.txt" "./LCL_x_up_1kb/" "lcl_exp_auto_genes.txt"
submit_job "lcl_x_resp_down_genes_auto.txt" "./LCL_x_down_1kb/" "lcl_exp_auto_genes.txt"
submit_job "lcl_y_resp_up_genes_auto.txt" "./LCL_y_up_1kb/" "lcl_exp_auto_genes.txt"
submit_job "lcl_y_resp_down_genes_auto.txt" "./LCL_y_down_1kb/" "lcl_exp_auto_genes.txt"
submit_job "lcl_x_y_resp_up_genes_auto.txt" "./LCL_x_y_up_1kb/" "lcl_exp_auto_genes.txt"
submit_job "lcl_x_y_resp_down_genes_auto.txt" "./LCL_x_y_down_1kb/" "lcl_exp_auto_genes.txt"

submit_job "fib_x_resp_up_genes_auto.txt" "./Fib_x_up_1kb/" "fib_exp_auto_genes.txt"
submit_job "fib_x_resp_down_genes_auto.txt" "./Fib_x_down_1kb/" "fib_exp_auto_genes.txt"
submit_job "fib_y_resp_up_genes_auto.txt" "./Fib_y_up_1kb/" "fib_exp_auto_genes.txt"
submit_job "fib_y_resp_down_genes_auto.txt" "./Fib_y_down_1kb/" "fib_exp_auto_genes.txt"
submit_job "fib_x_y_resp_up_genes_auto.txt" "./Fib_x_y_up_1kb/" "fib_exp_auto_genes.txt"
submit_job "fib_x_y_resp_down_genes_auto.txt" "./Fib_x_y_down_1kb/" "fib_exp_auto_genes.txt"
