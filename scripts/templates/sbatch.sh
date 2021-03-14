#!/bin/bash
#SBATCH --nodes=1                       # number of nodes requested
#SBATCH --ntasks=80                     # number of tasks (default: 1)
#SBATCH --partition=all                 # partition to run in (all or maxwell)
#SBATCH --job-name=NAME                 # job name
#SBATCH --output=NAME-%N-%j.out         # output file name
#SBATCH --error=NAME-%N-%j.err          # error file name
#SBATCH --time=4:00:00                  # runtime requested
#SBATCH --mail-user=ayan.paul@desy.de   # notification email
#SBATCH --mail-type=END,FAIL            # notification type
export LD_PRELOAD=""


# install the virtual environment and activate it
python3 -m venv .env
source .env/bin/activate

# install the necessary packages
pip3 install --upgrade pip
pip3 install -r requirements.txt

# run the classifier
python3 classifier.py

# deactivate and delete the environment
deactivate
rm -rf .env
