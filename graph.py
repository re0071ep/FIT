import numpy as np
import scipy.optimize
import matplotlib.pylab as plt
import math
from scipy import integrate
import cmath
from scipy.special import iv
import csv
from glob import glob
import os


def FitResult(freq, Z, Zfit_tot, chipName, file_i):
    plt.scatter(freq, Z, label="experiment",color="blue")
    plt.loglog(freq, Zfit_tot, label="fitting",color="red")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("|Z| [Ω]")
    plt.legend(loc='upper right', fontsize=10)
    plt.savefig("Result_" + str(chipName) + "/" + "FittingResult" + str(file_i))
    plt.close()

def EachFitResult_wCell(freq, freq_high, freq_mid, freq_low, Z, Z_high, Z_mid, Z_low, chipName, file_i):
    plt.scatter(freq, Z, label="experiment",color="blue")
    plt.loglog(freq_high, Z_high, label="fitting",color="red")
    plt.loglog(freq_mid, Z_mid, label="fitting",color="red")
    plt.loglog(freq_low, Z_low, label="fitting",color="red")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("|Z| [Ω]")
    plt.legend(loc='upper right', fontsize=10)
    plt.savefig("Result_" + str(chipName) + "/" + "HighMidLow Fitting" + str(file_i))
    plt.close()

def EachFitResult_woCell(freq, freq_high, freq_low, Z, Z_high, Z_low, chipName, file_i):
    plt.scatter(freq, Z, label="experiment",color="blue")
    plt.loglog(freq_high, Z_high, label="fitting",color="red")
    plt.loglog(freq_low, Z_low, label="fitting",color="red")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("|Z| [Ω]")
    plt.legend(loc='upper right', fontsize=10)
    plt.savefig("Result_" + str(chipName) + "/" + "HighMidLow Fitting" + str(file_i))
    plt.close()

def ErrorRsult(freq, Z, Zfit_tot, chipName, file_i):
    Z0 = freq/freq - 1
    Z_error = Zfit_tot*100 / Z - 100
    plt.scatter(freq, Z_error, label="Error rate",color="red")
    plt.plot(freq, Z0, label="",color="black")
    plt.xscale('log')
    plt.xlabel("Frequency[Hz]")
    plt.ylabel("Error rate [%]")
    plt.legend(loc='lower left', fontsize=10)
    plt.savefig("Result_" + str(chipName) + "/" + "error" + str(file_i))
    plt.close()

def Cdl0_result(exp_num, Cdl0_fit, chipName):
    plt.scatter(exp_num, Cdl0_fit, label="Cdl0",color="red")
    plt.xlabel("Experiment number")
    plt.ylabel("Cdl0 [F/m^2]")
    plt.legend(loc='upper right', fontsize=10)
    plt.savefig("Result_" + str(chipName) + "/" + "01_Fitting_Cdl0")
    plt.close()

def cellNum_result(exp_num, cellNum_fit, chipName):
    plt.scatter(exp_num, cellNum_fit, label="Cell number",color="red")
    plt.xlabel("Experiment number")
    plt.ylabel("Cell number")
    plt.legend(loc='upper right', fontsize=10)
    plt.savefig("Result_" + str(chipName) + "/" + "03_Fitting_CellNum")
    plt.close()

def cellRadius_result(exp_num, radius_fit, chipName):
    plt.scatter(exp_num, radius_fit*1e6, label="Cell radius",color="red")
    plt.xlabel("Experiment number")
    plt.ylabel("Cell radius [μm]")
    plt.legend(loc='upper right', fontsize=10)
    plt.savefig("Result_" + str(chipName) + "/" + "04_Fitting_CellRad")
    plt.close()

def Occupancy_result(exp_num, occupancy_fit, chipName):
    plt.scatter(exp_num, occupancy_fit, label="Occupancy",color="red")
    plt.xlabel("Experiment number")
    plt.ylabel("Occupancy [%]")
    plt.legend(loc='upper right', fontsize=10)
    plt.savefig("Result_" + str(chipName) + "/" + "05_Fitting_Occupancy")
    plt.close()

def R_result_wCell(exp_num, R_fit, chipName):
    plt.scatter(exp_num, R_fit, label="Rsol+Rgap+Rpar",color="red")
    plt.xlabel("Experiment number")
    plt.ylabel("Rsol+Rgap+Rpar [Ω]")
    plt.legend(loc='upper right', fontsize=10)
    plt.savefig("Result_" + str(chipName) + "/" + "02_Fitting_Rsol+Rgap+Rpar")
    plt.close()

def R_result_woCell(exp_num, R_fit, chipName):
    plt.scatter(exp_num, R_fit, label="Rsol+Rgap+Rpar",color="red")
    plt.xlabel("Experiment number")
    plt.ylabel("Rsol+Rpar [Ω]")
    plt.legend(loc='upper right', fontsize=10)
    plt.savefig("Result_" + str(chipName) + "/" + "02_Fitting_Rsol+Rpar")
    plt.close()