#%%
import numpy as np
import scipy.optimize
import matplotlib.pylab as plt
from scipy.special import iv
import csv
from glob import glob
import os
import graph
import FuncList

"""理論値"""
A = 2e-5 #電極面積
sigmasol = 1.88 #siraganian buffer の導電率　at 37℃
r_cell = 7e-6 #細胞半径
file_i = 0 #file数
file_num = 0 #file数

deviceName = input("Input device name (ALS = 0, LCRmeter = 1):\n")
Name = int(deviceName)
if Name == 0:
    device = "ALS"
elif Name == 1:
    device = "LCR"
Name = str(deviceName)
FuncList.kimura(Name)
date = input("Input experiment date (Example: 2019/11/27 → 191127): \n")
sampleName = input("Sample name: \n")
chipName = input("Chip name: \n")
fileType = input("File type (csv or txt): \n")
skipRows = input("Number of skip rows: \n")
usecols_freq = input("\nInput colum number of frequency, impedance and phase.\n\nColumn number of frequency: ")
usecols_imp = input("Column number of impedance: ")
usecols_phase = input("Column number of phase: ")
print("\n")
deviceName, usecols_freq, usecols_imp, usecols_phase, skipRows = int(deviceName), int(usecols_freq) - 1, int(usecols_imp) - 1, int(usecols_phase) - 1, int(skipRows)
low_range, high_range = 5, 10

path = os.getcwd() #このpyファイルのディレクトリを取得
os.mkdir(".\\Result_" + str(chipName))

txt = "Device name: " + device + "\nDate: " + str(date) + "\nSample　Name: " + str(sampleName) + "\nChip name: " + str(chipName) + "\nFile type (csv or txt): " + str(fileType) + "\nSkip row: " + str(skipRows) + "\nColumn number of frequency: " + str(usecols_freq) + "\nColumn number of impedance: " + str(usecols_imp) + "\nColumn number of phase: " + str(usecols_phase)
with open("Result_" + str(chipName) + "/00_InputParameter.txt", 'w') as t:
    t.write(txt)

with open("Result_" + str(chipName) + "/03_FittingResult.csv", 'w', newline="") as f: #fitting結果を出力するcsvfileの作成
    writer = csv.writer(f)
    writer.writerow(["num", "Cdl0[F/m^2]", "R", "Lpar"])


    
for file in glob(path + "/*" + str(date) + "**" + str(sampleName) + "*" + "_" + str(chipName) + "." + str(fileType)): #指定したフォルダ内のtxtファイルの読み込み
    file_num = file_num + 1
for file in glob(path + "/*" + str(date) + "*" + str(sampleName) + "_" + str(chipName) + "." + str(fileType)): #指定したフォルダ内のtxtファイルの読み込み
    file_i = file_i + 1
    data_result = 100 #消さないでください（最小２乗誤差のデータを抽出する際に使います）

    f = open(str(file))
    data = np.loadtxt(file, delimiter=',', usecols=[usecols_freq, usecols_imp, usecols_phase],skiprows = skipRows) #txtファイル内のデータを読み込み
    freq, Z, phase = np.array(data[:,0]), np.array(data[:,1]), np.array(data[:,2])

    para_oder = FuncList.paraOder(freq, Z, phase, deviceName)
    freq, Z, phase = para_oder[0], para_oder[1], para_oder[2]

    """Fittng frequency range"""
    phase_num = sum(1 for line in phase)
    FreqRange = FuncList.FreqAreaLowHigh_woCell(freq, Z, phase, low_range, high_range, phase_num)
    freq_low, freq_high, Z_low, Z_high = FreqRange[0], FreqRange[1], FreqRange[2], FreqRange[3]
    """Fitting low frequency"""
    Low_fit = FuncList.fittingLowFreq(A, freq_low, Z_low)
    Cdl0, p, Zfit_low = Low_fit[0], Low_fit[1], Low_fit[2]
    """fitting high frequency"""
    for high_sweep in range(2):
        high_fin = high_range + high_sweep
        high_start = 0 + high_sweep
        freq_high = np.array(freq[high_start : high_fin])
        Z_high = np.array(Z[high_start : high_fin])
        High_fit = FuncList.fittingHighFreq(A, freq_high, Z_high)
        R, Lpar, Zfit_high = High_fit[0], High_fit[1], High_fit[2]
        High_fit = FuncList.fitResult_woCell(freq, Z, A, Cdl0, p, Lpar, R, data_result)
        error_ave, data_result, Zfit_tot, phase_fit = High_fit[0], High_fit[1], High_fit[2], High_fit[3]

    print("Already finished: ", file_i, "/", file_num)
    with open("Result_" + str(chipName) + "/03_FittingResult.csv", 'a', newline="") as f:
        writer = csv.writer(f)
        writer.writerow([round(int(file_i),0), round(Cdl0,4), round(R,2), Lpar])
    graph.FitResult(freq, Z, Zfit_tot, chipName, file_i)
    graph.EachFitResult_woCell(freq, freq_high, freq_low, Z, Z_high, Z_low, chipName, file_i)
    graph.ErrorRsult(freq, Z, Zfit_tot, chipName, file_i)

    with open("Result_" + str(chipName) + "/gnuplot" + str(file_i) + ".csv", 'w', newline="") as gnu: #fitting結果を出力するcsvfileの作成
        writer = csv.writer(gnu)
        writer.writerow(["freq", "Z_exp[Ω]", "phase[-degree]", "Z_fit[Ω]", "phase_fit[-degree]"])
    with open("Result_" + str(chipName) + "/gnuplot" + str(file_i) + ".csv", 'a', newline="") as gnu:
        for i in range(phase_num):
            writer = csv.writer(gnu)
            writer.writerow([freq[i], round(Z[i], 3), round(phase[i], 3), round(Zfit_tot[i], 3), round(phase_fit[i], 3)])

fit_data = np.loadtxt("Result_"  + str(chipName) + "/03_FittingResult.csv", delimiter=',', skiprows = 1) #txtファイル内のデータを読み込み
print("\nExp.Num Cdl0[F/m^2]  Rsol+Rgap[Ω]  Lpar[H]\n", fit_data)
exp_num, Cdl0_fit, R_fit, Lpar_fit = np.array(fit_data[:,0]), np.array(fit_data[:,1]), np.array(fit_data[:,2]), np.array(fit_data[:,3])

graph.Cdl0_result(exp_num, Cdl0_fit, chipName)
graph.R_result_woCell(exp_num, R_fit, chipName)

input("\nPush Enter key\n")