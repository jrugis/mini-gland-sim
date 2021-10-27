import os
import time
import subprocess
import sys
import re

##################################################################
# functions
##################################################################

# make a directory-naming-friendly label from a parameter string
def make_label(s):
  tokens = s.strip().replace(".", "p").split()
  return "-".join(tokens)

##################################################################
# main program
##################################################################

print("mini_gland_sim")
run_dir = os.getcwd()

if(len(sys.argv) < 4):
  print("error: missing argument(s)")
  print("usage: python run_sim.py <mesh-zip-file> <parameter-file> <slurm-script>")
  quit()

mesh_zip_file = sys.argv[1] # mesh file
if not os.path.exists(run_dir + "/" + mesh_zip_file):
  print("no such mesh zip file: " + run_dir + "/" + mesh_zip_file)
  quit()

parms = sys.argv[2] # parameters file
if not os.path.exists(run_dir + "/" + parms):
  print("no such parameter file: " + parms)
  quit()

slurm = sys.argv[3] # Slurm script
if not os.path.exists(run_dir + "/" + slurm):
  print("no such Slurm script: " + slurm)
  quit()

# create the top level results directory
run_dir = re.sub("^/scale_wlg_persistent/filesets", "/nesi", run_dir)
dlist = run_dir.split("/")[-4:]
assert "project" in run_dir, "Run directory should be on project filesystem"
results_dir = run_dir.replace("project", "nobackup", 1)
results_dir = results_dir.replace("run", "results", 1)
results_dir += "/" + time.strftime("%y%m%d_%H%M%S")
print("results dir:", results_dir)
os.system("mkdir -p " + results_dir)
os.chdir(results_dir)

# copy some files into results directory
os.system("cp " + run_dir + "/mini-gland-sim .")
os.system("chmod 770 mini-gland-sim")
os.system("mkdir meshes")
os.system("unzip -q " + run_dir + "/" + mesh_zip_file + " -d ./meshes")
os.system("cp " + run_dir + "/" + parms + " p.ini")
os.system("cp " + run_dir + "/plot_*.py .")
os.system("cp " + run_dir + "/" + slurm + " run.sl")

# run the executable
#os.system("./mini-gland-sim")
os.system("sbatch run.sl")

# go back to top level
os.chdir(run_dir)
