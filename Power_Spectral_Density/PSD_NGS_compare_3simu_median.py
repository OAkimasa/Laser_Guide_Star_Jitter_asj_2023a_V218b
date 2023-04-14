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
    ngs_data = [row for row in reader]
    ngs_data = np.array(ngs_data, dtype=np.float64)  # obs
# t[ms], xslope_ave[pix], yslope_ave[pix], tip_coeff[arb], tilt_coeff[arb], focus_coeff[arb]
#print("ngs_data.shape:", ngs_data.shape)  # (48123, 6)

## シミュレーションデータの読み込み(mat)
ngs_simu_slopes_ave_median = scipy.io.loadmat('slopes_ave_update20230305_median.mat')  # slopes_ave median
ngs_simu_slopes_ave_median = ngs_simu_slopes_ave_median['slopes_ave']
ngs_simu_xslope_ave_median = ngs_simu_slopes_ave_median[:, 0]  # median
ngs_simu_yslope_ave_median = ngs_simu_slopes_ave_median[:, 1]  # median
ngs_simu_slopes_ave1wfs_median = scipy.io.loadmat('slopes_ave1wfs_update20230305_median.mat')  # slopes_ave median
ngs_simu_slopes_ave1wfs_median = ngs_simu_slopes_ave1wfs_median['slopes_ave1wfs']
ngs_simu_xslope_ave1wfs_median = ngs_simu_slopes_ave1wfs_median[:, 0]  # median
ngs_simu_yslope_ave1wfs_median = ngs_simu_slopes_ave1wfs_median[:, 1]  # median
ngs_simu_wave_tt_median = scipy.io.loadmat('wave_tt_update20230305_median.mat')  # slopes_ave median
ngs_simu_wave_tt_median = ngs_simu_wave_tt_median['wave_tt']
ngs_simu_wave_tip_median = ngs_simu_wave_tt_median[:, 0]  # median
ngs_simu_wave_tilt_median = ngs_simu_wave_tt_median[:, 1]  # median



## 観測データの切り出し
ngs_time = ngs_data[:, 0]*1e-3
ngs_simu_time = np.arange(1/400, (len(ngs_simu_xslope_ave_median)+1)/400, 1/400)  # range(400Hz*time)

## 観測データの処理
ngs_obs_xslope_ave = ngs_data[:, 1]  # [pix]
ngs_obs_yslope_ave = ngs_data[:, 2]  # [pix]
ngs_obs_xslope_ave = ngs_obs_xslope_ave - np.mean(ngs_obs_xslope_ave)  # obs [pix]
ngs_obs_yslope_ave = ngs_obs_yslope_ave - np.mean(ngs_obs_yslope_ave)  # obs [pix]

## シミュレーションデータの処理
ngs_simu_xslope_ave_median = ngs_simu_xslope_ave_median - np.mean(ngs_simu_xslope_ave_median)  # median [pix]
ngs_simu_yslope_ave_median = ngs_simu_yslope_ave_median - np.mean(ngs_simu_yslope_ave_median)  # median [pix]
ngs_simu_xslope_ave1wfs_median = ngs_simu_xslope_ave1wfs_median - np.mean(ngs_simu_xslope_ave1wfs_median)  # median [arcsec]
ngs_simu_yslope_ave1wfs_median = ngs_simu_yslope_ave1wfs_median - np.mean(ngs_simu_yslope_ave1wfs_median)  # median [arcsec]
ngs_simu_wave_tip_median = ngs_simu_wave_tip_median - np.mean(ngs_simu_wave_tip_median)  # median [arcsec]
ngs_simu_wave_tilt_median = ngs_simu_wave_tilt_median - np.mean(ngs_simu_wave_tilt_median)  # median [arcsec]


## unit: [pix] -> [arcsec]
ngs_obs_xslope_ave = ngs_obs_xslope_ave*slope_unit_conversion  # obs [arcsec]
ngs_obs_yslope_ave = ngs_obs_yslope_ave*slope_unit_conversion  # obs [arcsec]

ngs_simu_xslope_ave_median = ngs_simu_xslope_ave_median*simu_unit_conversion  # median [arcsec]
ngs_simu_yslope_ave_median = ngs_simu_yslope_ave_median*simu_unit_conversion  # median [arcsec]


## パワースペクトル密度の計算
fs = 400  # サンプリング周波数[Hz]

pw_ngs_obs_xslope = signal.welch(ngs_obs_xslope_ave, fs)  # obs
pw_ngs_obs_yslope = signal.welch(ngs_obs_yslope_ave, fs)  # obs

pw_ngs_simu_xslope_median = signal.welch(ngs_simu_xslope_ave_median, fs)  # median
pw_ngs_simu_yslope_median = signal.welch(ngs_simu_yslope_ave_median, fs)  # median
pw_ngs_simu_xslope1wfs_median = signal.welch(ngs_simu_xslope_ave1wfs_median, fs)  # median
pw_ngs_simu_yslope1wfs_median = signal.welch(ngs_simu_yslope_ave1wfs_median, fs)  # median
pw_ngs_simu_wave_tip_median = signal.welch(ngs_simu_wave_tip_median, fs)  # median
pw_ngs_simu_wave_tilt_median = signal.welch(ngs_simu_wave_tilt_median, fs)  # median

## G-tilt, Z-tiltの計算  Glenn A. Tyler 1992
f = np.linspace(0, 400, 400)  # [Hz]
G_tilt = (1e-3)*(f)**(-11/3)
Z_tilt = (1e-2)*(f)**(-17/3)


## データのプロット
fig_1 = plt.figure(figsize=(9, 9))
ax1 = plt.subplot(3, 1, 1)
ax1.plot(ngs_time, ngs_obs_xslope_ave, color="k", linewidth=0.2, label="NGS_obs_xslope")
ax1.plot(ngs_simu_time, ngs_simu_xslope_ave_median, color="green", linewidth=0.2, label="NGS_simu_xslope_median")
ax1.plot(ngs_simu_time, ngs_simu_xslope_ave1wfs_median, color="red", linewidth=0.2, label="NGS_simu_xslope1wfs_median")
ax1.plot(ngs_simu_time, ngs_simu_wave_tip_median, color="blue", linewidth=0.2, label="NGS_simu_wave_tip_median")
ax1.set_ylabel("xslope [arcsec]")
ax1.set_xlabel("Time [s]")
ax1.set_xlim(0, 10)
ax1.legend()

ax2 = plt.subplot(3, 1, 2)
ax2.plot(ngs_time, ngs_obs_yslope_ave, color="k", linewidth=0.2, label="NGS_obs_tilt")
ax2.plot(ngs_simu_time, ngs_simu_yslope_ave_median, color="green", linewidth=0.2, label="NGS_simu_tilt_median")
ax2.plot(ngs_simu_time, ngs_simu_yslope_ave1wfs_median, color="red", linewidth=0.2, label="NGS_simu_1_lensArray_tilt_median")
ax2.plot(ngs_simu_time, ngs_simu_wave_tilt_median, color="blue", linewidth=0.2, label="NGS_zernike_tilt_median")
ax2.set_ylabel("yslope [arcsec]")
ax2.set_xlabel("Time [s]")
ax2.set_xlim(0, 10)
ax2.legend()

# パワースペクトル密度のプロット
ax3 = plt.subplot(3, 1, 3)
ax3.plot(pw_ngs_obs_xslope[0], pw_ngs_obs_xslope[1], color="k", linewidth=1, label="NGS_obs_tip/tilt")
ax3.plot(pw_ngs_obs_yslope[0], pw_ngs_obs_yslope[1], color="k", linewidth=1, linestyle="dashed")
ax3.plot(pw_ngs_simu_xslope_median[0], pw_ngs_simu_xslope_median[1], color="green", linewidth=1, label="NGS_simu_tip/tilt_median")
ax3.plot(pw_ngs_simu_yslope_median[0], pw_ngs_simu_yslope_median[1], color="green", linewidth=1, linestyle="dashed")
ax3.plot(pw_ngs_simu_xslope1wfs_median[0], pw_ngs_simu_xslope1wfs_median[1], color="red", linewidth=1, label="NGS_simu_1_lensArray_tip/tilt_median")
ax3.plot(pw_ngs_simu_yslope1wfs_median[0], pw_ngs_simu_yslope1wfs_median[1], color="red", linewidth=1, linestyle="dashed")
ax3.plot(pw_ngs_simu_wave_tip_median[0], pw_ngs_simu_wave_tip_median[1], color="blue", linewidth=1, label="NGS_simu_zernike_tip/tilt_median")
ax3.plot(pw_ngs_simu_wave_tilt_median[0], pw_ngs_simu_wave_tilt_median[1], color="blue", linewidth=1, linestyle="dashed")
ax3.plot(f, G_tilt, color="orange", linewidth=1, label="G_tilt")
ax3.plot(f, Z_tilt, color="orange", linewidth=1, linestyle="dashdot", label="Z_tilt")
ax3.set_xlabel("Frequency [Hz]")
ax3.set_ylabel("PSD [arcsec^2/Hz]")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlim(1, 200)
#ax3.set_ylim(1e-6, 1e+2)
ax3.legend()

fig_1 = plt.figure(figsize=(9, 3.5))
# パワースペクトル密度のプロット
ax3 = plt.subplot(111)
ax3.plot(pw_ngs_obs_xslope[0], pw_ngs_obs_xslope[1], color="k", linewidth=1, label="NGS_obs_tip/tilt")
ax3.plot(pw_ngs_obs_yslope[0], pw_ngs_obs_yslope[1], color="k", linewidth=1, linestyle="dashed")
ax3.plot(pw_ngs_simu_xslope_median[0], pw_ngs_simu_xslope_median[1], color="green", linewidth=1, label="NGS_simu_tip/tilt_median")
ax3.plot(pw_ngs_simu_yslope_median[0], pw_ngs_simu_yslope_median[1], color="green", linewidth=1, linestyle="dashed")
ax3.plot(pw_ngs_simu_xslope1wfs_median[0], pw_ngs_simu_xslope1wfs_median[1], color="red", linewidth=1, label="NGS_simu_1_lensArray_tip/tilt_median")
ax3.plot(pw_ngs_simu_yslope1wfs_median[0], pw_ngs_simu_yslope1wfs_median[1], color="red", linewidth=1, linestyle="dashed")
ax3.plot(pw_ngs_simu_wave_tip_median[0], pw_ngs_simu_wave_tip_median[1], color="blue", linewidth=1, label="NGS_simu_before_lensArray_tip/tilt_median")
ax3.plot(pw_ngs_simu_wave_tilt_median[0], pw_ngs_simu_wave_tilt_median[1], color="blue", linewidth=1, linestyle="dashed")
ax3.plot(f, G_tilt, color="orange", linewidth=1, label="G_tilt")
ax3.plot(f, Z_tilt, color="orange", linewidth=1, linestyle="dashdot", label="Z_tilt")
ax3.set_xlabel("Frequency [Hz]")
ax3.set_ylabel("PSD [arcsec^2/Hz]")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlim(1, 200)
#ax3.set_ylim(1e-6, 1e+2)
#ax3.set_xlim(1e-1, 1e+2)
#ax3.set_ylim(1e-6, 1e-3)
ax3.legend()

plt.show()