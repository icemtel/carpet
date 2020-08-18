# Run jobs for different carpet sizes
# conda activate py37

cd "nx=10_ny=10"
python master_evals.py 2 1e-2 1e-8

cd "../nx=18_ny=18"
python master_evals.py 2 1e-2 1e-8
