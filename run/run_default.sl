#!/bin/bash -e
#SBATCH --account=nesi00119
#SBATCH --job-name=minigland
#SBATCH --time=01:00:00        # may need to be adjusted 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4      # could be adjusted
#SBATCH --hint=nomultithread
#SBATCH --mem=1G               # may need to be adjusted
#SBATCH --export=NONE          # these two lines required to make sure we get
export SLURM_EXPORT_ENV=ALL    # a clean environment in our Slurm job

# setup environment modules
ml purge
ml SUNDIALS/5.7.0-gimkl-2020a
ml h5pp/1.9.0-gimpi-2020a
ml Python/3.8.2-gimkl-2020a

# OpenMP configuration
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_NESTED=true
echo "Num threads: ${OMP_NUM_THREADS}"

# run the simulation
./mini-gland-sim

# plot some results
python plot_dev.py
python plot_electroneutrality.py
