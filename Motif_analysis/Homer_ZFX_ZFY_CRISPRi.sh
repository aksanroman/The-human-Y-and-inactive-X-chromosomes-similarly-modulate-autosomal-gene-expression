#!/bin/bash

# A script to run homer for all permutations.

# Job submission function
submit_job() {
    local input_file="$1"
    local output_dir="$2"


    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=homer
#SBATCH --output=homer_%j.out.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --partition=20    
#SBATCH --mail-type=END              

findMotifs.pl "$input_file" human "$output_dir" -bg fib_exp_auto_genes.txt -p 40 -start -1000 -end 1000 -nogo -fdr 100
EOF
}

# Submit jobs
submit_job "XX_zfx_v_cntrl_up_auto_genes.txt" "./XX_zfx_cntrl_up/"
submit_job "XX_zfx_v_cntrl_down_auto_genes.txt" "./XX_zfx_cntrl_down/" 
submit_job "XY_zfx_v_cntrl_up_auto_genes.txt" "./XY_zfx_cntrl_up/"
submit_job "XY_zfx_v_cntrl_down_auto_genes.txt" "./XY_zfx_cntrl_down/"
submit_job "XY_zfy_v_cntrl_up_auto_genes.txt" "./XY_zfy_cntrl_up/" 
submit_job "XY_zfy_v_cntrl_down_auto_genes.txt" "./XY_zfy_cntrl_down/"
submit_job "XY_zfxzfy_v_cntrl_up_auto_genes.txt" "./XY_zfxzfy_cntrl_up/"
submit_job "XY_zfxzfy_v_cntrl_down_auto_genes.txt" "./XY_zfxzfy_cntrl_down/"
