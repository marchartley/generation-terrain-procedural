import time
from multiprocessing import Process
import os

def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

def f(name):
    # info('function f')
    # print('hello', name, flush=False)
    time.sleep(3)
    print(f"Done {name}")

if __name__ == '__main__':
    # info('main line')
    start = time.time()
    ps = [Process(target=f, args=(i,)) for i in range(100)]
    for p in ps:
        p.start()
    for p in ps:
        p.join()
    print(time.time() - start)