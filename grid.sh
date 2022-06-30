#!/bin/bash
#SBATCH --mem-per-cpu=1024M
#SBATCH --time=03:00:00

julia --project=/home/aorfi/projects/def-raymond/aorfi/qHMC/qHMC_env ParameterSearch.jl