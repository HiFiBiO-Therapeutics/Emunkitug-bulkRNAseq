# Analysis_template
Default repo that contains commonly used files and structure.

## How to use this template

Create the repo, using this repo as a template. To take advantage of the structure of the repo, be sure
to run using a venv/conda environment with make and lintr enabled. Conda with Mamba installation is
recommended. Also recommended to have `git lfs` enabled

## Structure

### src

Folder that contains scripts that are used across different analysis notebooks

### data

Source data that the notebooks can call, including tables, csvs, rds, etc.

### output

Folder for analysis output. Usually used for graphs, images, and other things that won't be reladed for further analyses.

### makefile

Automated commands, including outputting conda environment, some git commands, and mapping of NAS servers if necessary

### .gitignore

Common files to ignore from git commits.

