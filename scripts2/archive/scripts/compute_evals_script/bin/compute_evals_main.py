# import multiprocessing as mp
import carpet.parallel_with_threads as pwt
import carpet
import logging
import subprocess

dirname = '.'
script_name = "compute_evals_script"
carpet.setup_logging('{}/evals.log'.format(dirname))


def main(k1,k2):
    logging.info("Starting: k1={} k2={}".format(k1, k2))
    completedProcess = subprocess.run(["python","-m",script_name, str(k1),str(k2)], cwd=dirname)
    if completedProcess.returncode == 1:
        raise RuntimeError("Subprocess finished with an error")
    logging.info("Finished: k1={} k2={}".format(k1, k2))

#, stdout=devnull, stderr=None
if __name__ == "__main__":
    nx = 6
    ny = 6
    list_of_args = [(k1,k2) for k1 in range(nx) for k2 in range(ny)]

    pwt.run_parallel(2, main, list_of_args)
