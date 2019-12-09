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
chipName = input("Chip name: \n")
fileType = input("File type (csv or txt): \n")
skipRows = input("Number of skip rows: \n")
usecols_freq = input("\nPlease input colum number of frequency, impedance and phase.\n\nColumn number of frequency: ")
usecols_imp = input("Column number of impedance: ")
usecols_phase = input("Column number of phase: ")
print("\n")
deviceName, usecols_freq, usecols_imp, usecols_phase, skipRows = int(deviceName), int(usecols_freq) - 1, int(usecols_imp) - 1, int(usecols_phase) - 1, int(skipRows)
i_low_range, i_high_range = 10, 10

path = os.getcwd() #このpyファイルのディレクトリを取得
os.mkdir(".\\Result_" + str(chipName))

txt = "Device name: " + device + "\nDate: " + str(date) +  "\nChip name: " + str(chipName) + "\nFile type (csv or txt): " + str(fileType) + "\nSkip row: " + str(skipRows) + "\nColumn number of frequency: " + str(usecols_freq) + "\nColumn number of impedance: " + str(usecols_imp) + "\nColumn number of phase: " + str(usecols_phase)
with open("Result_" + str(chipName) + "/00_InputParameter.txt", 'w') as t:
    t.write(txt)


with open("Result_" + str(chipName) + "/06_FittingResult.csv", 'w', newline="") as f: #fitting結果を出力するcsvfileの作成
    writer = csv.writer(f)
    writer.writerow(["num", "Cdl0[F/m^2]", "cell num", "cell radius[m]", "occupancy[%]", "R"])
for file in glob(path + "/*" + str(date) + "*" + "*_" + str(chipName) + "." + str(fileType)): 
    file_num = file_num + 1
for file in glob(path + "/*" + str(date) + "*" + "*_" + str(chipName) + "." + str(fileType)): 
    data_result = 100 #消さないでください
    file_i = file_i + 1
    f = open(str(file))
    data = np.loadtxt(file, delimiter=',', usecols=[usecols_freq, usecols_imp, usecols_phase],skiprows = skipRows) #txtファイル内のデータを読み込み
    freq, Z, phase = np.array(data[:,0]), np.array(data[:,1]), np.array(data[:,2])
    para_oder = FuncList.paraOder(freq, Z, phase, deviceName)
    freq, Z, phase = np.array(data[:,0]), np.array(data[:,1]), np.array(data[:,2])
    """Fittng frequency range"""
    phase_num = sum(1 for line in phase)
    FreqRange = FuncList.FreqAreaLowHigh_wCell(freq, Z, phase, i_low_range, i_high_range, phase_num)
    freq_low, freq_high, Z_low, Z_high = FreqRange[0], FreqRange[1], FreqRange[2], FreqRange[3]
    """Fitting low frequency"""
    Low_fit = FuncList.fittingLowFreq(A, freq_low, Z_low)
    Cdl0, p, Zfit_low = Low_fit[0], Low_fit[1], Low_fit[2]
    """fitting high frequency"""
    High_fit = FuncList.fittingHighFreq(A, freq_high, Z_high)
    R, Lpar, Zfit_high = High_fit[0], High_fit[1], High_fit[2]
    """fitting middle frequency"""
    for i_mid_range in range(10):
        i_mid_range = i_mid_range + 2
        for mid_sweep in range(20):
            FreqRange = FuncList.FreqAreaLowHigh_wCell(freq, Z, phase, i_low_range, i_high_range, phase_num)
            num_mid = FreqRange[4]
            midFreq_start = num_mid - i_mid_range - mid_sweep + 15
            midFreq_fin = num_mid + i_mid_range - mid_sweep + 15
            freq_mid = np.array(freq[midFreq_start : midFreq_fin])
            Z_mid = np.array(Z[midFreq_start : midFreq_fin])
            if file_i < 2: #１つめのデータ
                Mid_fit = FuncList.fittingMidFreq1(A, freq_mid, Z_mid, sigmasol, Lpar, r_cell, p, Cdl0, R)
                hgap, cell_num, Acell, Zfit_mid = Mid_fit[0], Mid_fit[1], Mid_fit[2], Mid_fit[3]
                fitResult = FuncList.fitResult1_wCell(freq, Z, A, cell_num, r_cell, Cdl0, p, sigmasol, hgap, Lpar, R, data_result)
                error_ave, data_result, Zfit_tot, phase_fit = fitResult[0], fitResult[1], fitResult[2], fitResult[3]
            else: #2つめ以降のデータ
                Mid_fit = FuncList.fittingMidFreq2(A, freq_mid, Z_mid, sigmasol, Lpar, cell_num, p, Cdl0, R)
                hgap, r_cell, Acell, Zfit_mid = Mid_fit[0], Mid_fit[1], Mid_fit[2], Mid_fit[3]
                fitResult = FuncList.fitResult2_wCell(freq, Z, A, cell_num, r_cell, Cdl0, p, sigmasol, hgap, Lpar, R, data_result)
                error_ave, data_result, Zfit_tot, phase_fit = fitResult[0], fitResult[1], fitResult[2], fitResult[3]
    print("Already finished: ", file_i, "/", file_num)
    with open("Result_" + str(chipName) + "/06_FittingResult.csv", 'a', newline="") as f:
        writer = csv.writer(f)
        writer.writerow([round(file_i,0), round(Cdl0, 3), round(cell_num, 0), round(r_cell, 9), round(100*Acell/A, 1), round(R, 2)])

    graph.FitResult(freq, Z, Zfit_tot, chipName, file_i)
    graph.EachFitResult_wCell(freq, freq_high, freq_mid, freq_low, Z, Z_high, Z_mid, Z_low, chipName, file_i)
    graph.ErrorRsult(freq, Z, Zfit_tot, chipName, file_i)

    with open("Result_" + str(chipName) + "/gnuplot" + str(file_i) + ".csv", 'w', newline="") as gnu: #fitting結果を出力するcsvfileの作成
        writer = csv.writer(gnu)
        writer.writerow(["freq", "Z_exp[Ω]", "phase[-degree]", "Z_fit[Ω]", "phase_fit[-degree]"])
    with open("Result_" + str(chipName) + "/gnuplot" + str(file_i) + ".csv", 'a', newline="") as gnu:
        for i in range(phase_num):
            writer = csv.writer(gnu)
            writer.writerow([freq[i], round(Z[i], 3), round(phase[i], 3), round(Zfit_tot[i], 3), round(phase_fit[i], 3)])

fit_data = np.loadtxt("Result_"  + str(chipName) + "/06_FittingResult.csv", delimiter=',', skiprows = 1) #txtファイル内のデータを読み込み
print("\nExp.Num Cdl0[F/m^2] CellNum[cells] CellRadius[m] Occupancy[%] Rsol+Rgap+Rpar\n", fit_data)
exp_num, Cdl0_fit, cellNum_fit, radius_fit, occupancy_fit, R_fit = np.array(fit_data[:,0]), np.array(fit_data[:,1]), np.array(fit_data[:,2]), np.array(fit_data[:,3]), np.array(fit_data[:,4]), np.array(fit_data[:,5])

graph.Cdl0_result(exp_num, Cdl0_fit, chipName)
graph.cellNum_result(exp_num, cellNum_fit, chipName)
graph.cellRadius_result(exp_num, radius_fit, chipName)
graph.Occupancy_result(exp_num, occupancy_fit, chipName)
graph.R_result_wCell(exp_num, R_fit, chipName)

input("\nPush Enter key\n")