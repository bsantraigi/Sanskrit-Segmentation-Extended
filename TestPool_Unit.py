from multiprocessing import Process
import multiprocessing as mp
import os
from sentences import *
import numpy as np
from Train_n_Save_NNet import *

def pooled_Test(modelFile, vpid, queue, testFolder, filePerProcess = 100 ):
    print('Child process with vpid:{}, pid:{} started.'.format(vpid, os.getpid()))
    trainer = Trainer()
    trainer.Load(modelFile)

    loaded_SKT = pickle.load(open('../Simultaneous_CompatSKT_10K.p', 'rb'))
    loaded_DCS = pickle.load(open('../Simultaneous_DCS_10K.p', 'rb'))

    TestFiles = []
    for f in os.listdir(testFolder):
        if '.ds.bz2' in f:
            TestFiles.append(f)

    loader = pickle.load(open('../bz2Dataset_10K.p', 'rb'))
    TestFiles = loader['TestFiles']
    TrainFiles = loader['TrainFiles']

    for i in range(vpid*filePerProcess, vpid*filePerProcess + filePerProcess):
        fn = TestFiles[i]
        fn = fn.replace('.ds.bz2', '.p2')

        dsbz2_name = testFolder + TestFiles[i]

        sentenceObj = loaded_SKT[fn]
        dcsObj = loaded_DCS[fn]
        results = trainer.Test(sentenceObj, dcsObj, dsbz2_name)
        if results is not None:
            queue.put(results)

    print('Child process with vpid:{}, pid:{} closed.'.format(vpid, os.getpid()))
