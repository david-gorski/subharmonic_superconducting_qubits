import numpy as np
import matplotlib.pyplot as plt
import h5py
import fittingAndDataProcess as fdp
from IQdata import loadH5IntoIQData  #HaPiCodes is a package written by ourselves
import pickle

ssbArray = np.linspace(0.085, 0.175, 91) #sideband frequency we used to modulate the qubit generator's frequency
timeArray = np.linspace(30, 230, 51, dtype=int)  #length of the smoothbox pulse

saveDir = rf"/Users/david/Code/subharmonic_superconducting_qubits/mingkang_graphing/TimeRabi_amp0.15_time30to230_51_SSB0.085to0.175_91/"  #file saved directory

IList = np.zeros((len(ssbArray), len(timeArray)))
QList = np.zeros((len(ssbArray), len(timeArray)))

if __name__ == '__main__':
    for i, iSSB in enumerate(ssbArray):  # read all the data from files. the files format are HDF5, it can be processed with python package h5py

        fileLoad = loadH5IntoIQData(saveDir, f'TimeRabi_time_ssb{i}_{iSSB}')
        IQdata = fileLoad[0]
        paramDict = fileLoad[1]

        ssbArray = paramDict['ssbArray']
        timeArray = paramDict['timeArray']

        Id, Qd = IQdata.I_rot, IQdata.Q_rot
        I_avg, Q_avg = fdp.average_data(Id, Qd)
        iqNew = fdp.rotate_complex(I_avg, Q_avg, 0)
        I_avg = iqNew.real
        Q_avg = iqNew.imag

        IList[i] = I_avg
        QList[i] = Q_avg
    iqNew = fdp.rotate_complex(IList, QList, 0.9)
    I_rot, Q_rot = iqNew.real, iqNew.imag
    plt.figure(figsize=(14, 5))
    plt.subplot(1,2,1)
    plt.pcolormesh(1.4195098 - ssbArray, timeArray, I_rot.T, shading='auto')
    plt.xlabel("freq(GHz)")
    plt.ylabel("pulse length(ns)")
    plt.colorbar()
    plt.subplot(1,2,2)
    plt.pcolormesh(1.4195098 - ssbArray, timeArray, Q_rot.T, shading='auto')
    plt.xlabel("freq(GHz)")
    plt.ylabel("pulse length(ns)")
    plt.colorbar()
    plt.show()

    plt.figure(figsize=(14, 5))
    plt.subplot(1,2,1)
    plt.pcolormesh(ssbArray, timeArray, I_rot.T, shading='auto')
    plt.xlabel("SSB freq(GHz)")
    plt.ylabel("pulse length(ns)")
    plt.colorbar()
    plt.subplot(1,2,2)
    plt.pcolormesh(ssbArray, timeArray, Q_rot.T, shading='auto')
    plt.xlabel("SSB freq(GHz)")
    plt.ylabel("pulse length(ns)")
    plt.colorbar()
    plt.show()

    name = "mingkang_data_0.15" + ".pkl"
    data = {
        "freqArray": 1.4195098 - ssbArray,
        "ssbArray": ssbArray,
        "timeArray": timeArray,
        "Q_rot.T": Q_rot.T,
        "I_rot.T": I_rot.T
    }
    with open(name, "wb") as f:
        pickle.dump(data, f)