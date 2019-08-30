#import multiprocessing as mp
import carpet.parallel_with_threads as pwt
import carpet
import logging
import subprocess

dirname = '.'
script_name = "compute_evals_script"
carpet.setup_logging('{}/evals.log'.format(dirname))

# su

def main(k1,k2):
    logging.info("Starting: k1={} k2={}".format(k1, k2))
    completedProcess = subprocess.run(["python","-m",script_name, str(k1),str(k2)], cwd=dirname)
    logging.info("Finished: k1={} k2={}".format(k1, k2))

#, stdout=devnull, stderr=None
if __name__ == "__main__":
    nx = 6
    ny = 6
    list_of_args = [(0,0),(1,1)]#[(k1,k2) for k1 in range(nx) for k2 in range(ny)]

    pwt.run_parallel(2, main, list_of_args)

    # def calc_sort_evals_evecs(L_log_lin):
    #     # Ljapunov spectrum
    #     evals, evecs = sp.linalg.eig(L_log_lin)
    #
    #     # Sort eigenvalues
    #     idx = evals.argsort()
    #     evals = evals[idx]
    #     evecs = evecs[:, idx]
    #     return evals, evecs
    #
    #
    # def calc_evals_evecs_mtwist(k1, k2, delta0, tol):
    #     '''
    #     evecs: eigenvectors as columns!
    #     '''
    #     L = get_L_mtwist(k1, k2, delta0, tol)
    #     L_log_lin = L - sp.eye(N)
    #     return calc_sort_evals_evecs(L_log_lin)
