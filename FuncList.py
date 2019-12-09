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

"""Freq, Imp, Phaseの列指定&列の順番決め"""
def paraOder(freq, Z, phase, deviceName):
    if deviceName == 1:
        """LCRでの取得データはALSデータの行を反転させた形なので、行の反転を行う"""
        freq = np.flipud(freq)
        Z = np.flipud(Z)
        phase = np.flipud(phase)
    """低周波部のfreq, imp, phaseの大小関係を用い、列の順番を決める"""
    row_num = sum(1 for line in freq)
    for para_order in range(5):
        if np.array(freq[row_num-1]) > np.array(Z[row_num-1]):
            freq_data = Z
            Z_data = freq
            Z = Z_data
            freq = freq_data
        elif np.array(phase[row_num-1]) > np.array(Z[row_num-1]):
            Z_data = phase
            phase_data = Z
            Z = Z_data
            phase = phase_data
        elif np.array(phase[row_num-1]) > np.array(freq[row_num-1]):
            freq_data = phase
            phase_data = freq
            phase = phase_data
            freq = freq_data
    return(freq, Z, phase)

"""Fittingに用いる周波数範囲の指定（細胞あり）"""
def FreqAreaLowHigh_wCell(freq, Z, phase, i_low_range, i_high_range, phase_num):
    """low frequency"""
    num = sum(1 for line in freq) - 2
    lowFreq_start = num - i_low_range
    lowFreq_fin = num
    freq_low = np.array(freq[lowFreq_start : lowFreq_fin])
    Z_low = np.array(Z[lowFreq_start : lowFreq_fin])
    """high frequency"""
    highFreq_fin = i_high_range
    highFreq_start = 0
    freq_high = np.array(freq[highFreq_start : highFreq_fin])
    Z_high =np.array(Z[highFreq_start : highFreq_fin])
    for n_phase in range(phase_num):
        phase_base = np.array(phase[n_phase]) + 45
        if phase_base > 0:
            if phase_base < 10:
                num_mid = n_phase + 1
                break
    return(freq_low, freq_high, Z_low, Z_high, num_mid)

"""Fittingに用いる周波数範囲の指定（細胞なし）"""
def FreqAreaLowHigh_woCell(freq, Z, phase, i_low_range, i_high_range, phase_num):
    """Fitting low frequency"""
    num = sum(1 for line in freq) - 2
    lowFreq_start = num - i_low_range
    lowFreq_fin = num
    freq_low = np.array(freq[lowFreq_start : lowFreq_fin])
    Z_low = np.array(Z[lowFreq_start : lowFreq_fin])
    """fitting high frequency"""
    highFreq_fin = i_high_range #fitting範囲の左端指定
    highFreq_start = 0 #fitting範囲の右端指定
    freq_high = np.array(freq[highFreq_start : highFreq_fin])
    Z_high =np.array(Z[highFreq_start : highFreq_fin])
    return(freq_low, freq_high, Z_low, Z_high)

"""Cdl0, pのパラメータ抽出"""
def fittingLowFreq(A, freq_low, Z_low):
    parameter_initial1 = np.array([0.1, 0.5]) #Cdl0, pの初期値指定
    parameter_bounds1 = ([0,0],[1,1]) #フィッティング値の範囲
    def func1(freq_low, Cdl0, p):
        Zdl  = 1 / (((2j*np.pi*freq_low)**p)*Cdl0*A)
        Ztot = 2*Zdl
        return abs(Ztot)
    paramater_optimal1, covariance = scipy.optimize.curve_fit(func1, freq_low, Z_low, p0=parameter_initial1, bounds=parameter_bounds1)
    fitting = func1(Z_low,paramater_optimal1[0],paramater_optimal1[1])
    Cdl0 = paramater_optimal1[0]
    p = paramater_optimal1[1]
    Zfit_low = abs(2 / (((2j*np.pi*freq_low)**p)*Cdl0*A))
    return(Cdl0, p, Zfit_low)

"""Rsol+Rgap+Rparのパラメータ抽出"""
def fittingHighFreq(A, freq_high, Z_high):
    parameter_initial2 = np.array([20, 2e-5]) #Rsol+Rgap+Rpar, Lparの初期値指定
    parameter_bounds2 = ([0,0],[100, 1])#フィッティング値の範囲
    def func2(freq_high, R, Lpar):
        Zl = 2j*np.pi*freq_high*Lpar
        Ztot =  R +  Zl
        return abs(Ztot)
    paramater_optimal2, covariance = scipy.optimize.curve_fit(func2, freq_high, Z_high, p0=parameter_initial2, bounds=parameter_bounds2)
    fitting = func2(Z_high,paramater_optimal2[0],paramater_optimal2[1])
    R = paramater_optimal2[0]
    Lpar = paramater_optimal2[1]
    Zfit_high = abs(R + 2j*np.pi*freq_high*Lpar)
    return(R, Lpar, Zfit_high)

"""隙間インピのフィッティング（抗原抗体反応前の１回目の測定データ使用）"""
def fittingMidFreq1(A, freq_mid, Z_mid, sigmasol, Lpar, r_cell, p, Cdl0, R):
    parameter_initial3 = np.array([50e-9, 10000]) #hgap, cell_numの初期値指定
    parameter_bounds3 = ([0,1],[100e-9,1000000000]) #フィッティング値の範囲
    def func3(freq_mid, hgap, cell_num):
        Zdl0 = 1 / (((2j*np.pi*freq_mid)**p)*Cdl0)
        c = ((1/(sigmasol*hgap)) * (1/Zdl0))**(1/2)
        I0 = iv(0, c*r_cell)
        I1 = iv(1, c*r_cell)
        Acell = cell_num*np.pi*(r_cell)**2
        Zc = (Zdl0/(Acell)*c*r_cell*I0) / (2*I1)
        Zl = 2j*np.pi*freq_mid*Lpar
        Ztot = 2*(Zc*(Zdl0/(A-Acell)))/(Zc+Zdl0/(A-Acell)) + R + Zl
        return abs(Ztot)
    paramater_optimal3, covariance = scipy.optimize.curve_fit(func3, freq_mid, Z_mid, p0=parameter_initial3, bounds=parameter_bounds3, maxfev = 5000000)
    fitting = func3(Z_mid,paramater_optimal3[0],paramater_optimal3[1])
    hgap = paramater_optimal3[0]
    cell_num = paramater_optimal3[1]
    Acell = cell_num*np.pi*(r_cell)**2
    Zdl_ele  = 1 / (((2j*np.pi*freq_mid)**p)*Cdl0*(A-Acell))
    Zdl_cell = 1 / (((2j*np.pi*freq_mid)**p)*Cdl0*Acell)
    Zdl0 = 1 / (((2j*np.pi*freq_mid)**p)*Cdl0)
    c = ((1/(sigmasol*hgap)) * (1/Zdl0))**(1/2)
    I0 = iv(0, c*r_cell)
    I1 = iv(1, c*r_cell)
    Zc = (Zdl0/(Acell)*c*r_cell*I0) / (2*I1)
    Zl = 2j*np.pi*freq_mid*Lpar
    Zfit_mid = abs(2*(Zc*(Zdl0/(A-Acell)))/(Zc+Zdl0/(A-Acell)) + R + Zl)
    return(hgap, cell_num, Acell, Zfit_mid)

"""隙間インピのフィッティング（抗原抗体反応前の１回目以降の測定データ使用）"""
def fittingMidFreq2(A, freq_mid, Z_mid, sigmasol, Lpar, cell_num, p, Cdl0, R):
    parameter_initial4 = np.array([50e-9, 7e-6]) #hgap, r_cellの初期値指定
    parameter_bounds4 = ([0,0],[100e-9,1]) #フィッティング値の範囲
    def func4(freq_mid, hgap, r_cell_fit):
        Zdl0 = 1 / (((2j*np.pi*freq_mid)**p)*Cdl0)
        c = ((1/(sigmasol*hgap)) * (1/Zdl0))**(1/2)
        I0 = iv(0, c*r_cell_fit)
        I1 = iv(1, c*r_cell_fit)
        Acell = cell_num*np.pi*(r_cell_fit)**2
        Zc = (Zdl0/(Acell)*c*r_cell_fit*I0) / (2*I1)
        Zl = 2j*np.pi*freq_mid*Lpar
        Ztot = 2*(Zc*(Zdl0/(A-Acell)))/(Zc+Zdl0/(A-Acell)) + R + Zl
        return abs(Ztot)
    paramater_optimal4, covariance = scipy.optimize.curve_fit(func4, freq_mid,   Z_mid, p0=parameter_initial4, bounds=parameter_bounds4, maxfev = 5000000)
    fitting = func4(Z_mid,paramater_optimal4[0],paramater_optimal4[1])
    hgap = paramater_optimal4[0]
    r_cell = paramater_optimal4[1]
    Acell = cell_num*np.pi*(r_cell)**2
    Zdl_ele  = 1 / (((2j*np.pi*freq_mid)**p)*Cdl0*(A-Acell))
    Zdl_cell = 1 / (((2j*np.pi*freq_mid)**p)*Cdl0*Acell)
    Zdl0 = 1 / (((2j*np.pi*freq_mid)**p)*Cdl0)
    c = ((1/(sigmasol*hgap)) * (1/Zdl0))**(1/2)
    I0 = iv(0, c*r_cell)
    I1 = iv(1, c*r_cell)
    Zc = (Zdl0/(Acell)*c*r_cell*I0) / (2*I1)
    Zl = 2j*np.pi*freq_mid*Lpar
    Zfit_mid = abs(2*(Zc*(Zdl0/(A-Acell)))/(Zc+Zdl0/(A-Acell)) + R)
    return(hgap, r_cell, Acell, Zfit_mid)

"""パラメータ抽出結果を用いてfittingインピーダンス&実験値との誤差を計算"""
def fitResult_woCell(freq, Z, A, Cdl0, p, Lpar, R, data_result):
    Zfit_tot = 2 / (((2j*np.pi*freq)**p)*Cdl0*A) + R + 2j*np.pi*freq*Lpar
    Zfit_re = Zfit_tot.real
    Zfit_im = Zfit_tot.imag
    phase_fit = np.rad2deg(np.arctan(Zfit_im / Zfit_re))
    Zfit_tot = abs(2 / (((2j*np.pi*freq)**p)*Cdl0*A) + R + 2j*np.pi*freq*Lpar)
    error_ave = sum(((Zfit_tot/Z)-1)**2) / sum(1 for line in freq)
    fit_data = round(error_ave, 4)
    """2乗平均誤差の一番少ないfitting結果を選択"""
    if fit_data < data_result:
        data_result = fit_data
    return(error_ave, data_result, Zfit_tot, phase_fit)

"""パラメータ抽出結果を用いてfittingインピーダンス&実験値との誤差を計算"""
"""抗原抗体反応前の１回目の測定データ使用"""
def fitResult1_wCell(freq, Z, A, cell_num, r_cell, Cdl0, p, sigmasol, hgap, Lpar, R, data_result):
    Acell = cell_num*np.pi*(r_cell)**2
    Zdl_ele  = 1 / (((2j*np.pi*freq)**p)*Cdl0*(A-Acell))
    Zdl_cell = 1 / (((2j*np.pi*freq)**p)*Cdl0*Acell)
    Zdl0 = 1 / (((2j*np.pi*freq)**p)*Cdl0)
    c = ((1/(sigmasol*hgap)) * (1/Zdl0))**(1/2)
    I0 = iv(0, c*r_cell)
    I1 = iv(1, c*r_cell)
    Zc = (Zdl0/(Acell)*c*r_cell*I0) / (2*I1)
    Zl = 2j*np.pi*freq*Lpar
    Zfit_tot = 2*(Zc*(Zdl0/(A-Acell)))/(Zc+Zdl0/(A-Acell)) + R + Zl
    Zfit_re = Zfit_tot.real
    Zfit_im = Zfit_tot.imag
    Zfit_tot = abs(2*(Zc*(Zdl0/(A-Acell)))/(Zc+Zdl0/(A-Acell)) + R + Zl)
    phase_fit = np.rad2deg(np.arctan(Zfit_im / Zfit_re))
    error_ave = sum(((Zfit_tot/Z)-1)**2) / sum(1 for line in freq)
    fit_data = round(error_ave, 4)
    """2乗平均誤差の一番少ないfitting結果を選択"""
    if fit_data < data_result:
        data_result = fit_data
    return(error_ave, data_result, Zfit_tot, phase_fit)

"""パラメータ抽出結果を用いてfittingインピーダンス&実験値との誤差を計算"""
"""抗原抗体反応前の２回目以降の測定データ使用"""
def fitResult2_wCell(freq, Z, A, cell_num, r_cell, Cdl0, p, sigmasol, hgap, Lpar, R, data_result):
    Acell = cell_num*np.pi*(r_cell)**2
    Zdl_ele  = 1 / (((2j*np.pi*freq)**p)*Cdl0*(A-Acell))
    Zdl_cell = 1 / (((2j*np.pi*freq)**p)*Cdl0*Acell)
    Zdl0 = 1 / (((2j*np.pi*freq)**p)*Cdl0)
    c = ((1/(sigmasol*hgap)) * (1/Zdl0))**(1/2)
    I0 = iv(0, c*r_cell)
    I1 = iv(1, c*r_cell)
    Zc = (Zdl0/(Acell)*c*r_cell*I0) / (2*I1)
    Zl = 2j*np.pi*freq*Lpar
    Zfit_tot = 2*(Zc*(Zdl0/(A-Acell)))/(Zc+Zdl0/(A-Acell)) + R + Zl
    Zfit_re = Zfit_tot.real
    Zfit_im = Zfit_tot.imag
    Zfit_tot = abs(2*(Zc*(Zdl0/(A-Acell)))/(Zc+Zdl0/(A-Acell)) + R + Zl)
    phase_fit = np.rad2deg(np.arctan(Zfit_im / Zfit_re))
    error_ave = sum(((Zfit_tot/Z)-1)**2) / sum(1 for line in freq)
    fit_data = round(error_ave, 4)
    if fit_data < data_result:
        data_result = fit_data
    return(error_ave, data_result, Zfit_tot, phase_fit)

def kimura(Name):
    if Name == "0107":
        print("\n----------------------------------------------------------------")
        print("This fitting tool was produced by KIMURA Kaiken at 2019/11/28")
        print("----------------------------------------------------------------\n")
        print("2019/11/28     Ver. 1.0 ")
        print("2019/12/09     Ver. 1.1 ")
        print("")

