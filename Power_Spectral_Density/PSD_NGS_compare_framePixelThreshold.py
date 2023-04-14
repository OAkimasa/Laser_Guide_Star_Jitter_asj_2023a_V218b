import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.signal as signal
import scipy.io


# 1.0[pix]=0.57[arcsec], tip_coeff(tilt_coeff) 109=4.0[pix]=2.28[arcsec]
slope_unit_conversion = 0.57  # arcsec/pix
coeff_unit_conversion = 2.28/109.  # ~0.0209 arcsec/pix
simu_unit_conversion = ((589.*1e-9)/(8./(32.)))*(1./(2.*2.))*(180./np.pi)*3600.  # ~0.121 arcsec/pix
print("simu_unit_conversion:", simu_unit_conversion)


## 観測データの読み込み(csv)
with open('20221110_ngs_open_400Hz_40_TTFlist.csv') as f:
    reader = csv.reader(f)
    ngs_obs = [row for row in reader]
    ngs_obs = np.array(ngs_obs, dtype=np.float64)
# t[ms], xslope_ave[pix], yslope_ave[pix], tip_coeff[arb], tilt_coeff[arb], focus_coeff[arb]
#print("ngs_obs.shape:", ngs_obs.shape)  # (48123, 6)


## シミュレーションデータの読み込み(mat)
ngs_simu_slopes_ave = scipy.io.loadmat('slopes_ave_update20230303.mat')
ngs_simu_slopes_ave = ngs_simu_slopes_ave['slopes_ave']
ngs_simu_xslope_ave = ngs_simu_slopes_ave[:, 0]
ngs_simu_yslope_ave = ngs_simu_slopes_ave[:, 1]

ngs_simu_slopes_ave_threshold3e4 = scipy.io.loadmat('slopes_ave_update20230303_threshold3e4.mat')
ngs_simu_slopes_ave_threshold3e4 = ngs_simu_slopes_ave_threshold3e4['slopes_ave']
ngs_simu_xslope_ave_threshold3e4 = ngs_simu_slopes_ave_threshold3e4[:, 0]
ngs_simu_yslope_ave_threshold3e4 = ngs_simu_slopes_ave_threshold3e4[:, 1]

ngs_simu_slopes_ave_threshold4e4 = scipy.io.loadmat('slopes_ave_update20230303_threshold4e4.mat')
ngs_simu_slopes_ave_threshold4e4 = ngs_simu_slopes_ave_threshold4e4['slopes_ave']
ngs_simu_xslope_ave_threshold4e4 = ngs_simu_slopes_ave_threshold4e4[:, 0]
ngs_simu_yslope_ave_threshold4e4 = ngs_simu_slopes_ave_threshold4e4[:, 1]

ngs_simu_slopes_ave_threshold5e4 = scipy.io.loadmat('slopes_ave_update20230303_threshold5e4.mat')
ngs_simu_slopes_ave_threshold5e4 = ngs_simu_slopes_ave_threshold5e4['slopes_ave']
ngs_simu_xslope_ave_threshold5e4 = ngs_simu_slopes_ave_threshold5e4[:, 0]
ngs_simu_yslope_ave_threshold5e4 = ngs_simu_slopes_ave_threshold5e4[:, 1]

ngs_simu_slopes_ave_threshold6e4 = scipy.io.loadmat('slopes_ave_update20230303_threshold6e4.mat')
ngs_simu_slopes_ave_threshold6e4 = ngs_simu_slopes_ave_threshold6e4['slopes_ave']
ngs_simu_xslope_ave_threshold6e4 = ngs_simu_slopes_ave_threshold6e4[:, 0]
ngs_simu_yslope_ave_threshold6e4 = ngs_simu_slopes_ave_threshold6e4[:, 1]


## データの切り出し
ngs_time = ngs_obs[:, 0]*1e-3  # [ms] -> [s]
ngs_simu_time = np.arange(1/400, (len(ngs_simu_xslope_ave)+1)/400, 1/400)  # range(400Hz*time)


## データの処理
ngs_obs_xslope_ave = ngs_obs[:, 1]  # [pix]
ngs_obs_yslope_ave = ngs_obs[:, 2]  # [pix]
ngs_obs_xslope_ave = ngs_obs_xslope_ave - np.mean(ngs_obs_xslope_ave)  # [pix]
ngs_obs_yslope_ave = ngs_obs_yslope_ave - np.mean(ngs_obs_yslope_ave)  # [pix]

ngs_simu_xslope_ave = ngs_simu_xslope_ave - np.mean(ngs_simu_xslope_ave)  # OOMAO_framePixelThreshold0
ngs_simu_yslope_ave = ngs_simu_yslope_ave - np.mean(ngs_simu_yslope_ave)  # framePixelThreshold0

ngs_simu_xslope_ave_threshold3e4 = ngs_simu_xslope_ave_threshold3e4 - np.mean(ngs_simu_xslope_ave_threshold3e4)  # framePixelThreshold3e4
ngs_simu_yslope_ave_threshold3e4 = ngs_simu_yslope_ave_threshold3e4 - np.mean(ngs_simu_yslope_ave_threshold3e4)  # framePixelThreshold3e4

ngs_simu_xslope_ave_threshold4e4 = ngs_simu_xslope_ave_threshold4e4 - np.mean(ngs_simu_xslope_ave_threshold4e4)  # framePixelThreshold4e4
ngs_simu_yslope_ave_threshold4e4 = ngs_simu_yslope_ave_threshold4e4 - np.mean(ngs_simu_yslope_ave_threshold4e4)  # framePixelThreshold4e4

ngs_simu_xslope_ave_threshold5e4 = ngs_simu_xslope_ave_threshold5e4 - np.mean(ngs_simu_xslope_ave_threshold5e4)  # framePixelThreshold5e4
ngs_simu_yslope_ave_threshold5e4 = ngs_simu_yslope_ave_threshold5e4 - np.mean(ngs_simu_yslope_ave_threshold5e4)  # framePixelThreshold5e4

ngs_simu_xslope_ave_threshold6e4 = ngs_simu_xslope_ave_threshold6e4 - np.mean(ngs_simu_xslope_ave_threshold6e4)  # framePixelThreshold6e4
ngs_simu_yslope_ave_threshold6e4 = ngs_simu_yslope_ave_threshold6e4 - np.mean(ngs_simu_yslope_ave_threshold6e4)  # framePixelThreshold6e4


## unit: [pix] -> [arcsec]
ngs_obs_xslope_ave = ngs_obs_xslope_ave*slope_unit_conversion
ngs_obs_yslope_ave = ngs_obs_yslope_ave*slope_unit_conversion

ngs_simu_xslope_ave = ngs_simu_xslope_ave*simu_unit_conversion
ngs_simu_yslope_ave = ngs_simu_yslope_ave*simu_unit_conversion

ngs_simu_xslope_ave_threshold3e4 = ngs_simu_xslope_ave_threshold3e4*simu_unit_conversion
ngs_simu_yslope_ave_threshold3e4 = ngs_simu_yslope_ave_threshold3e4*simu_unit_conversion

ngs_simu_xslope_ave_threshold4e4 = ngs_simu_xslope_ave_threshold4e4*simu_unit_conversion
ngs_simu_yslope_ave_threshold4e4 = ngs_simu_yslope_ave_threshold4e4*simu_unit_conversion

ngs_simu_xslope_ave_threshold5e4 = ngs_simu_xslope_ave_threshold5e4*simu_unit_conversion
ngs_simu_yslope_ave_threshold5e4 = ngs_simu_yslope_ave_threshold5e4*simu_unit_conversion

ngs_simu_xslope_ave_threshold6e4 = ngs_simu_xslope_ave_threshold6e4*simu_unit_conversion
ngs_simu_yslope_ave_threshold6e4 = ngs_simu_yslope_ave_threshold6e4*simu_unit_conversion


## PSDの計算
fs = 400  # サンプリング周波数[Hz]
pw_ngs_obs_xslope = signal.welch(ngs_obs_xslope_ave, fs)
pw_ngs_obs_yslope = signal.welch(ngs_obs_yslope_ave, fs)

pw_ngs_simu_xslope = signal.welch(ngs_simu_xslope_ave, fs)
pw_ngs_simu_yslope = signal.welch(ngs_simu_yslope_ave, fs)

pw_ngs_simu_xslope_threshold3e4 = signal.welch(ngs_simu_xslope_ave_threshold3e4, fs)
pw_ngs_simu_yslope_threshold3e4 = signal.welch(ngs_simu_yslope_ave_threshold3e4, fs)

pw_ngs_simu_xslope_threshold4e4 = signal.welch(ngs_simu_xslope_ave_threshold4e4, fs)
pw_ngs_simu_yslope_threshold4e4 = signal.welch(ngs_simu_yslope_ave_threshold4e4, fs)

pw_ngs_simu_xslope_threshold5e4 = signal.welch(ngs_simu_xslope_ave_threshold5e4, fs)
pw_ngs_simu_yslope_threshold5e4 = signal.welch(ngs_simu_yslope_ave_threshold5e4, fs)

pw_ngs_simu_xslope_threshold6e4 = signal.welch(ngs_simu_xslope_ave_threshold6e4, fs)
pw_ngs_simu_yslope_threshold6e4 = signal.welch(ngs_simu_yslope_ave_threshold6e4, fs)


## データのプロット
fig_1 = plt.figure(figsize=(9, 9))
ax1 = plt.subplot(3, 1, 1)
ax1.plot(ngs_time, ngs_obs_xslope_ave, color="red", linewidth=0.2, label="NGS_obs_xslope")
ax1.plot(ngs_simu_time, ngs_simu_xslope_ave, color="green", linewidth=0.2, label="NGS_simu_xslope_framePixelThreshold0")
ax1.plot(ngs_simu_time[:4000], ngs_simu_xslope_ave_threshold3e4, color="blue", linewidth=0.2, label="NGS_simu_xslope_framePixelThreshold3e4")
ax1.plot(ngs_simu_time[:4000], ngs_simu_xslope_ave_threshold4e4, color="orange", linewidth=0.2, label="NGS_simu_xslope_framePixelThreshold4e4")
ax1.plot(ngs_simu_time[:4000], ngs_simu_xslope_ave_threshold5e4, color="purple", linewidth=0.2, label="NGS_simu_xslope_framePixelThreshold5e4")
ax1.plot(ngs_simu_time[:4000], ngs_simu_xslope_ave_threshold6e4, color="black", linewidth=0.2, label="NGS_simu_xslope_framePixelThreshold6e4")
ax1.set_ylabel("xslope [arcsec]")
ax1.set_xlabel("Time [s]")
ax1.set_xlim(0, 10)
ax1.legend()

ax2 = plt.subplot(3, 1, 2)
ax2.plot(ngs_time, ngs_obs_yslope_ave, color="red", linewidth=0.2, label="NGS_obs_yslope")
ax2.plot(ngs_simu_time, ngs_simu_yslope_ave, color="green", linewidth=0.2, label="NGS_simu_yslope_framePixelThreshold0")
ax2.plot(ngs_simu_time[:4000], ngs_simu_yslope_ave_threshold3e4, color="blue", linewidth=0.2, label="NGS_simu_yslope_framePixelThreshold3e4")
ax2.plot(ngs_simu_time[:4000], ngs_simu_yslope_ave_threshold4e4, color="orange", linewidth=0.2, label="NGS_simu_yslope_framePixelThreshold4e4")
ax2.plot(ngs_simu_time[:4000], ngs_simu_yslope_ave_threshold5e4, color="purple", linewidth=0.2, label="NGS_simu_yslope_framePixelThreshold5e4")
ax2.plot(ngs_simu_time[:4000], ngs_simu_yslope_ave_threshold6e4, color="black", linewidth=0.2, label="NGS_simu_yslope_framePixelThreshold6e4")
ax2.set_ylabel("yslope [arcsec]")
ax2.set_xlabel("Time [s]")
ax2.set_xlim(0, 10)
ax2.legend()

# パワースペクトル密度のプロット
ax3 = plt.subplot(3, 1, 3)
ax3.plot(pw_ngs_obs_xslope[0], pw_ngs_obs_xslope[1], color="red", linewidth=1, label="NGS_obs_xyslope")
ax3.plot(pw_ngs_obs_yslope[0], pw_ngs_obs_yslope[1], color="red", linewidth=1, linestyle="dashed")
ax3.plot(pw_ngs_simu_xslope[0], pw_ngs_simu_xslope[1], color="green", linewidth=1, label="NGS_simu_xyslope_framePixelThreshold0")
ax3.plot(pw_ngs_simu_yslope[0], pw_ngs_simu_yslope[1], color="green", linewidth=1, linestyle="dashed")
ax3.plot(pw_ngs_simu_xslope_threshold3e4[0], pw_ngs_simu_xslope_threshold3e4[1], color="blue", linewidth=1, label="NGS_simu_xyslope_framePixelThreshold3e4")
ax3.plot(pw_ngs_simu_yslope_threshold3e4[0], pw_ngs_simu_yslope_threshold3e4[1], color="blue", linewidth=1, linestyle="dashed")
ax3.plot(pw_ngs_simu_xslope_threshold4e4[0], pw_ngs_simu_xslope_threshold4e4[1], color="orange", linewidth=1, label="NGS_simu_xyslope_framePixelThreshold4e4")
ax3.plot(pw_ngs_simu_yslope_threshold4e4[0], pw_ngs_simu_yslope_threshold4e4[1], color="orange", linewidth=1, linestyle="dashed")
ax3.plot(pw_ngs_simu_xslope_threshold5e4[0], pw_ngs_simu_xslope_threshold5e4[1], color="purple", linewidth=1, label="NGS_simu_xyslope_framePixelThreshold5e4")
ax3.plot(pw_ngs_simu_yslope_threshold5e4[0], pw_ngs_simu_yslope_threshold5e4[1], color="purple", linewidth=1, linestyle="dashed")
ax3.plot(pw_ngs_simu_xslope_threshold6e4[0], pw_ngs_simu_xslope_threshold6e4[1], color="black", linewidth=1, label="NGS_simu_xyslope_framePixelThreshold6e4")
ax3.plot(pw_ngs_simu_yslope_threshold6e4[0], pw_ngs_simu_yslope_threshold6e4[1], color="black", linewidth=1, linestyle="dashed")
ax3.set_xlabel("Frequency [Hz]")
ax3.set_ylabel("PSD [arcsec^2/Hz]")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlim(1, 200)
#ax3.set_ylim(1e-6, 1e+2)
ax3.legend()

plt.show()