#!/bin/bash -e
#SBATCH --account=nesi00119
#SBATCH --job-name=minigland
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem=1G

ml purge
ml SUNDIALS/5.7.0-gimkl-2020a
ml h5pp/1.9.0-gimpi-2020a
ml Python/3.8.2-gimkl-2020a

./mini-gland-sim

python plot_dev.py
python plot_electroneutrality.py
