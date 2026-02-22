for todo in *.in; do
    # Define the name of the output SLURM file
    slm_file="${todo}.slm"
    
    # Create the SLURM script for each TODO value
    cat > "$slm_file" <<EOL
#!/bin/bash
#SBATCH --job-name=$todo
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --time=
#SBATCH --nodes=
#SBATCH -C cpu
#SBATCH --ntasks
#SBATCH --cpus-per-task
#SBATCH -e ${todo}.err
#SBATCH --qos=regular

module load qchem
export QCLOCALSCR=\$QCSCRATCH

###run q-chem
qchem  $todo ${todo%.in}.out

###

echo "Finished running $todo"
echo "================="
echo " "
echo "Job ended at :"
date
EOL

    #echo "Created SLURM file: $slm_file"
done

