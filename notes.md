# Mini-Gland Simulation on Mahuika

## Get the code

Only do these steps once when you first get started:

```sh
mkdir -p /nesi/project/nesi00119/$USER
cd /nesi/project/nesi00119/$USER
git clone https://github.com/jrugis/mini-gland-sim.git
cd mini-gland-sim
```

## Update the code to the latest version

Run the following to get the latest version of the code:

```sh
git pull
```

If the above command didn't work (maybe your local repo has diverged somehow),
you could run the following to reset to the upstream version:

```sh
git fetch origin
git reset --hard origin/master
```

## Prerequisites for building and running

*Note these could be put in ~/.bash_profile*

Load the required software environment modules.

```sh
ml purge
ml SUNDIALS/5.7.0-gimkl-2020a
ml h5pp/1.9.0-gimpi-2020a
ml Python/3.8.2-gimkl-2020a
```

## Build the code

Whenever changes are made to the C++ source code you must rebuild it.

```sh
cd mini-gland-sim/build
make clean
make -j
```

## Run the code

```sh
cd mini-gland-sim/run

cp parms_default.ini myparms.ini
# edit myparms.ini to set simulation parameters

cp run_default.sl myrun.sl
# edit Slurm script myrun.sl (e.g. update time limit, number of CPUs per task, etc.)

python run_sim_dev.py mini_gland_slurm-20373309_1.ply myparms.ini myrun.sl
# or more generally: python run_sim_dev.py <meshfile> <parameterfile> <slurmscript>

# check the Slurm queue
squeue -u $USER
```

## Results

Results can be found in the directory reported by the run script, it will look something like:

```
/nesi/nobackup/nesi00119/${USER}/mini-gland-sim/results/XXXX
```

You can follow the simulation progress by viewing the output files (*.out).
