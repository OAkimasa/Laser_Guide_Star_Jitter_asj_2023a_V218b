import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.signal as signal
import scipy.io


## unit conversion
# observation: 1.0[pix]=0.57[arcsec], tip_coeff(tilt_coeff) 109=4.0[pix]=2.28[arcsec]
slope_unit_conversion = 0.57  # arcsec/pix
# simulation: 4.0[m]=512[pix]
simu_unit_conversion = np.rad2deg(np.arcsin(1/(90*1000/np.cos(np.radians(90.0-55)))))*3600/128  # ~0.0147 arcsec/pix
#print("simu_unit_conversion:", simu_unit_conversion)


## 観測データの読み込み(csv)
with open('20221110_lgs_open_400Hz_3_TTFlist.csv') as f:  # LGS_observation
    reader = csv.reader(f)
    lgs_data = [row for row in reader]
    lgs_data_400Hz_EL55_obs = np.array(lgs_data, dtype=np.float64)
#print("lgs_data_400Hz_EL55_obs.shape", lgs_data_400Hz_EL55_obs.shape)  # (48121, 6)

with open('20221110_ngs_open_400Hz_40_TTFlist.csv') as f:  # NGS_observation
    reader = csv.reader(f)
    ngs_data = [row for row in reader]
    ngs_data = np.array(ngs_data, dtype=np.float64)
# t[ms], xslope_ave[pix], yslope_ave[pix], tip_coeff[arb], tilt_coeff[arb], focus_coeff[arb]
#print("ngs_data.shape:", ngs_data.shape)  # (48123, 6)


## シミュレーションデータの読み込み(csv)
with open('popt_list_200Hz_EL55_update_median.csv') as f:  # median
    reader = csv.reader(f)
    lgs_data = [row for row in reader]
    LGS_data_200Hz_EL55_simu_median = np.array(lgs_data, dtype=np.float64)  # median
#print("LGS_data_200Hz_EL55_simu.shape", LGS_data_200Hz_EL55_simu.shape)  # (6000, 6)
with open('popt_list_200Hz_EL55_update_good.csv') as f:  # good
    reader = csv.reader(f)
    lgs_data = [row for row in reader]
    LGS_data_200Hz_EL55_simu_good = np.array(lgs_data, dtype=np.float64)  # good
with open('popt_list_200Hz_EL55_update_bad.csv') as f:  # bad
    reader = csv.reader(f)
    lgs_data = [row for row in reader]
    LGS_data_200Hz_EL55_simu_bad = np.array(lgs_data, dtype=np.float64)  # bad


## lgs観測データ
# データの処理
lgs_time_obs = lgs_data_400Hz_EL55_obs[:, 0]*1e-3  # [s]
lgs_xslope_ave_obs = lgs_data_400Hz_EL55_obs[:, 1]  # [pix]
lgs_yslope_ave_obs = lgs_data_400Hz_EL55_obs[:, 2]  # [pix]
lgs_xslope_ave_obs = lgs_xslope_ave_obs - np.mean(lgs_xslope_ave_obs)  # [pix] obs
lgs_yslope_ave_obs = lgs_yslope_ave_obs - np.mean(lgs_yslope_ave_obs)  # [pix] obs
lgs_EL55_tip_obs = lgs_xslope_ave_obs*slope_unit_conversion  # [arcsec] obs
lgs_EL55_tilt_obs = lgs_yslope_ave_obs*slope_unit_conversion  # [arcsec] obs

# PSDの計算
fs = 400  # サンプリング周波数[Hz]
pw_lgs_EL55_tip_obs = signal.welch(lgs_EL55_tip_obs, fs) # obs
pw_lgs_EL55_tilt_obs = signal.welch(lgs_EL55_tilt_obs, fs) # obs


## ngs観測データ
# データの処理
ngs_time_obs = ngs_data[:, 0]*1e-3  # [s]
ngs_xslope_ave_obs = ngs_data[:, 1]  # [pix]
ngs_yslope_ave_obs = ngs_data[:, 2]  # [pix]
ngs_xslope_ave_obs = ngs_xslope_ave_obs - np.mean(ngs_xslope_ave_obs)  # [pix] obs
ngs_yslope_ave_obs = ngs_yslope_ave_obs - np.mean(ngs_yslope_ave_obs)  # [pix] obs
ngs_EL55_tip_obs = ngs_xslope_ave_obs*slope_unit_conversion  # [arcsec] obs
ngs_EL55_tilt_obs = ngs_yslope_ave_obs*slope_unit_conversion  # [arcsec] obs

# PSDの計算
fs = 400  # サンプリング周波数[Hz]
pw_ngs_EL55_tip_obs = signal.welch(ngs_EL55_tip_obs, fs) # obs
pw_ngs_EL55_tilt_obs = signal.welch(ngs_EL55_tilt_obs, fs) # obs


## lgsシミュレーション 200Hz注意
# データの処理
lgs_time_simu = np.arange(1/200, (len(LGS_data_200Hz_EL55_simu_median[:, 1])+1)/200, 1/200)  # [s], 200Hz, shape=(6000,)
lgs_EL55_x_pos_simu_median = LGS_data_200Hz_EL55_simu_median[:, 1]-256  # [pix] simu_median
lgs_EL55_y_pos_simu_median = LGS_data_200Hz_EL55_simu_median[:, 2]-256  # [pix] simu_median
lgs_EL55_tip_simu_median = lgs_EL55_x_pos_simu_median*simu_unit_conversion  # [arcsec] simu_median
lgs_EL55_tilt_simu_median = lgs_EL55_y_pos_simu_median*simu_unit_conversion  # [arcsec] simu_median
lgs_EL55_x_pos_simu_good = LGS_data_200Hz_EL55_simu_good[:, 1]-256  # [pix] simu_good
lgs_EL55_y_pos_simu_good = LGS_data_200Hz_EL55_simu_good[:, 2]-256  # [pix] simu_good
lgs_EL55_tip_simu_good = lgs_EL55_x_pos_simu_good*simu_unit_conversion  # [arcsec] simu_good
lgs_EL55_tilt_simu_good = lgs_EL55_y_pos_simu_good*simu_unit_conversion  # [arcsec] simu_good
lgs_EL55_x_pos_simu_bad = LGS_data_200Hz_EL55_simu_bad[:, 1]-256  # [pix] simu_bad
lgs_EL55_y_pos_simu_bad = LGS_data_200Hz_EL55_simu_bad[:, 2]-256  # [pix] simu_bad
lgs_EL55_tip_simu_bad = lgs_EL55_x_pos_simu_bad*simu_unit_conversion  # [arcsec] simu_bad
lgs_EL55_tilt_simu_bad = lgs_EL55_y_pos_simu_bad*simu_unit_conversion  # [arcsec] simu_bad


# PSDの計算
fs = 200  # サンプリング周波数[Hz]
pw_lgs_EL55_tip_simu_median = signal.welch(lgs_EL55_tip_simu_median, fs) # simu_median
pw_lgs_EL55_tilt_simu_median = signal.welch(lgs_EL55_tilt_simu_median, fs) # simu_median
pw_lgs_EL55_tip_simu_good = signal.welch(lgs_EL55_tip_simu_good, fs) # simu_good
pw_lgs_EL55_tilt_simu_good = signal.welch(lgs_EL55_tilt_simu_good, fs) # simu_good
pw_lgs_EL55_tip_simu_bad = signal.welch(lgs_EL55_tip_simu_bad, fs) # simu_bad
pw_lgs_EL55_tilt_simu_bad = signal.welch(lgs_EL55_tilt_simu_bad, fs) # simu_bad


## データのプロット
fig_1 = plt.figure(figsize=(9, 7))
"""ax1 = plt.subplot(3, 1, 1)
ax1.plot(lgs_time_obs, lgs_EL55_tip_obs, color="orange", linewidth=0.2, label="LGS_EL55_tip_obs")
ax1.plot(lgs_time_simu, lgs_EL55_tip_simu_bad, color="red", linewidth=0.2, label="LGS_EL55_tip_simu_bad")
ax1.plot(lgs_time_simu, lgs_EL55_tip_simu_median, color="blue", linewidth=0.2, label="LGS_EL55_tip_simu_median")
ax1.plot(lgs_time_simu, lgs_EL55_tip_simu_good, color="green", linewidth=0.2, label="LGS_EL55_tip_simu_good")
ax1.plot(ngs_time_obs, ngs_EL55_tip_obs, color="k", linewidth=0.2, label="NGS_EL55_tip_obs")
ax1.set_ylabel("Tip [arcsec]")
ax1.set_xlabel("Time [s]")
ax1.set_xlim(0, 10)
ax1.legend(loc='upper right', fontsize=10)"""

ax2 = plt.subplot(2, 1, 1)
ax2.plot(lgs_time_obs, lgs_EL55_tilt_obs, color="orange", linewidth=0.2, label="LGS_EL55_tilt_obs")
ax2.plot(lgs_time_simu, lgs_EL55_tilt_simu_bad, color="red", linewidth=0.2, label="LGS_EL55_tilt_simu_bad")
ax2.plot(lgs_time_simu, lgs_EL55_tilt_simu_median, color="blue", linewidth=0.2, label="LGS_EL55_tilt_simu_median")
ax2.plot(lgs_time_simu, lgs_EL55_tilt_simu_good, color="green", linewidth=0.2, label="LGS_EL55_tilt_simu_good")
ax2.plot(ngs_time_obs, ngs_EL55_tilt_obs, color="k", linewidth=0.2, label="NGS_EL55_tilt_obs")
ax2.set_ylabel("Tilt [arcsec]")
ax2.set_xlabel("Time [s]")
ax2.set_xlim(0, 10)
ax2.legend(loc='upper right', fontsize=10)

ax3 = plt.subplot(2, 1, 2)
ax3.plot(pw_lgs_EL55_tip_obs[0], pw_lgs_EL55_tip_obs[1], color="orange", linewidth=0.6, label="LGS_EL55_tip/tilt_obs")
ax3.plot(pw_lgs_EL55_tilt_obs[0], pw_lgs_EL55_tilt_obs[1], color="orange", linewidth=0.6, linestyle="dashed")
ax3.plot(pw_lgs_EL55_tip_simu_bad[0], pw_lgs_EL55_tip_simu_bad[1], color="red", linewidth=0.5, label="LGS_EL55_tip/tilt_simu_bad")
ax3.plot(pw_lgs_EL55_tilt_simu_bad[0], pw_lgs_EL55_tilt_simu_bad[1], color="red", linewidth=0.5, linestyle="dashed")
ax3.plot(pw_lgs_EL55_tip_simu_median[0], pw_lgs_EL55_tip_simu_median[1], color="blue", linewidth=0.5, label="LGS_EL55_tip/tilt_simu_median")
ax3.plot(pw_lgs_EL55_tilt_simu_median[0], pw_lgs_EL55_tilt_simu_median[1], color="blue", linewidth=0.5, linestyle="dashed")
ax3.plot(pw_lgs_EL55_tip_simu_good[0], pw_lgs_EL55_tip_simu_good[1], color="green", linewidth=0.5, label="LGS_EL55_tip/tilt_simu_good")
ax3.plot(pw_lgs_EL55_tilt_simu_good[0], pw_lgs_EL55_tilt_simu_good[1], color="green", linewidth=0.5, linestyle="dashed")
ax3.plot(pw_ngs_EL55_tip_obs[0], pw_ngs_EL55_tip_obs[1], color="k", linewidth=0.5, label="NGS_EL55_tip/tilt_obs")
ax3.plot(pw_ngs_EL55_tilt_obs[0], pw_ngs_EL55_tilt_obs[1], color="k", linewidth=0.5, linestyle="dashed")
ax3.set_xlabel("Frequency [Hz]")
ax3.set_ylabel("PSD [arcsec^2/Hz]")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlim(1, 200)
#ax3.set_ylim(1e-6, 1e+2)
ax3.legend(loc='lower left', fontsize=10)

plt.show()