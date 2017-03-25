from multiprocessing import Process
import multiprocessing as mp
import os
from sentences import *
import numpy as np

Trainer_Module = None
SKT_DCS_pairs = None

q = mp.Queue()
def InitModule(_trainer_module, _skt_dcs_pairs):
    global Trainer_Module, SKT_DCS_pairs
    Trainer_Module = _trainer_module
    SKT_DCS_pairs = _skt_dcs_pairs

def multi_wrapper(q, skt, dcs):
    TM = q.get()
    print('NPSUM: ', TM.Test(skt, dcs))

def inner_multiprocess():
    q.put(Trainer_Module)
    q.put(Trainer_Module)
    q.put(Trainer_Module)
    q.put(Trainer_Module)

    print(Trainer_Module)
    procs = [None]*4
    for i in range(4):
        procs[i] = Process(target = multi_wrapper, args = (q, SKT_DCS_pairs[i][0], SKT_DCS_pairs[i][1]))

    for proc in procs:
        proc.start()

    for proc in procs:
        proc.join()
