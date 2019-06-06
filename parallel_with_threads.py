'''
Copied from `hydro` package on 2019-06-06 - Anton Solovev.

Use `run_parallel` function to parallelize a job.
Before use make sure that the function does not change some objects which are shared by the other instances of the function.
'''

import threading
import queue

def report_progress_until_finished(q, time_between_reports, time_between_checks, qsize_init=None):
    import logging
    import time
    if qsize_init is None:
        qsize_init = q.qsize()
    last_report_time = time.time()
    while not q.empty():
        time_elapsed = time.time() - last_report_time
        if time_elapsed > time_between_reports:
            portion_finished = 1 - q.qsize() / qsize_init
            logging.info("run_parallel: {:.1%} jobs are finished".format(portion_finished))
            last_report_time = time.time()
        time.sleep(time_between_checks)

def run_parallel(num_processes, function_to_run, list_of_args=None, list_of_kwargs=None,
                 report_progress=True, time_between_reports=3600, time_between_checks=5):
    '''
    Function will be executed in parallel,
    taking each time arguments from list_of_args, and key-arguments from list_of_kwargs
    The functions exits only after each thread is finished.
    Limitations:
        (1) If there is one argument, it has to be wrapped in a list, e.g.
          [[1],[2],[3],..] instead of [1,2,3,4..]
        (2) Works with lists; not tested with generators?
    '''
    # Check that input is correct
    if list_of_args == None and list_of_kwargs == None:
        raise ValueError('Both args and kwargs are None')
    elif list_of_args == None:
        list_of_args = [[] for i in range(len(list_of_kwargs))]
    elif list_of_kwargs == None:
        list_of_kwargs = [{} for i in range(len(list_of_args))]
    elif len(list_of_args) != len(list_of_kwargs):
        raise ValueError('List of args and list of kwargs have different length')

    # Define a worker which will continuously execute function taking input parameters from the queue
    def worker():
        while True:
            next_argument = q.get()
            if next_argument is None:
                break
            args, kwargs = next_argument
            function_to_run(*args, **kwargs)
            q.task_done()

    # Initialize the queue
    q = queue.Queue()
    source = zip(list_of_args, list_of_kwargs)
    for item in source:
        q.put(item)
    qsize_init = q.qsize()
    # Initialize the threads
    threads = []
    for i in range(num_processes):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)
    # Report progress every hour if this option is enabled
    if report_progress:
        report_progress_until_finished(q, time_between_reports, time_between_checks, qsize_init)
    # Block until all tasks are done
    q.join()

    # stop workers
    for i in range(num_processes):
        q.put(None)
    for t in threads:
        t.join()


if __name__ == '__main__':
    import logging

    logging.basicConfig(level=20, format='{asctime} {levelname} {threadName} {name} {message}', style='{')

    # Testing the function
    a = iter(range(100))  # Test how common data is handled


    def f(x, y):
        import random
        import time
        t = 0.5 * random.randint(1, 10)
        print(threading.get_ident(), threading.current_thread().name)
        print(threading.current_thread(), "Next in a:", next(a), " Result:", x * y, 'sleep for', t)
        time.sleep(t)
        return x * x

    print(threading.get_ident(), threading.current_thread().name)

    list_of_args = None  # [[i, i ** 3] for i in range(7)]
    list_of_kwargs = [{'x': i, 'y': i * i} for i in range(20)]
    run_parallel(2, f, list_of_args, list_of_kwargs, time_between_reports=3, time_between_checks=1)  # [[k] for k in range(10)])
