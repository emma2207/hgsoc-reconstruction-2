#!/bin/bash
#SBATCH --job-name=experiment_nf
#SBATCH --account=amc-general
#SBATCH --output=experiment_nf_%J.log
#SBATCH --error=experiment_nf_%J.err
#SBATCH --time=08:00:00
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --mem=4G
#SBATCH --nodes=1 

set -eo pipefail

# --------------------------------------------------
# 0) Resolve project root (directory of this script)
# --------------------------------------------------
# Check if SLURM_SUBMIT_DIR is set (indicates Slurm environment)
if [ -n "${SLURM_SUBMIT_DIR}" ]; then
    # Running on HPC via Slurm
    PRJ_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    RUN_MODE="HPC"
    echo "Running on HPC (Slurm). Project directory: ${PRJ_DIR}"
else
    # Running locally
    # Get the directory where the script itself is located when run locally
    PRJ_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    RUN_MODE="LOCAL"
    echo "Running locally. Project directory: ${PRJ_DIR}"
fi

# Change to the project directory
cd "${PRJ_DIR}"
echo "Current working directory: $(pwd)"

# --------------------------------------------------
# 0.5) Specifying paths specific to Alpine
# --------------------------------------------------
# export paths specific to Alpine
if [ "${RUN_MODE}" == "HPC" ]; then
    export PATH=/usr/include:$PATH
    export CPATH=/usr/include/:$CPATH
    export C_INCLUDE_PATH=/usr/include/:$C_INCLUDE_PATH
fi

# --------------------------------------------------
# 1) Ensure that Nextflow is installed
# --------------------------------------------------
if [ "${RUN_MODE}" == "HPC" ]; then
    module load miniforge
    module load nextflow
fi

if ! command -v nextflow &> /dev/null
then
    echo "Error: Nextflow is not installed or not found in your PATH."
    echo "Please install Nextflow and ensure it's accessible in your system's PATH."
    echo "You can typically install it with: curl -s https://get.nextflow.io | bash"
    echo "Then move the 'nextflow' executable to a directory in your PATH (e.g., ~/bin or /usr/local/bin)."
    exit 1 # Exit the script with an error code
fi

# --------------------------------------------------
# 2) Run the pipeline
# --------------------------------------------------
echo "••• Launching Nextflow"

# Define Nextflow's work directory based on run mode
NEXTFLOW_WORK_DIR="${PRJ_DIR}/nextflow"

if [ "${RUN_MODE}" == "HPC" ]; then
    nextflow run main.nf -profile slurm -w "${NEXTFLOW_WORK_DIR}" -process.echo -resume -with-trace
else
    nextflow run main.nf -w "${NEXTFLOW_WORK_DIR}" -process.echo -resume -with-trace
fi

echo "••• Pipeline finished 🎉"
