#!/bin/bash
#SBATCH --mem-per-cpu=5G
#SBATCH --time=10:00:00

julia --project=/home/aorfi/projects/def-raymond/aorfi/qHMC/qHMC_env SKParameterSearch.jl