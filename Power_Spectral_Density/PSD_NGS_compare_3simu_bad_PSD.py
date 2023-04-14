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
    ngs_obs = np.array(ngs_obs, dtype=np.float64)  # obs
# t[ms], xslope_ave[pix], yslope_ave[pix], tip_coeff[arb], tilt_coeff[arb], focus_coeff[arb]
#print("ngs_obs.shape:", ngs_obs.shape)  # (48123, 6)

## シミュレーションデータの読み込み(mat)
ngs_simu_slopes_ave_bad = scipy.io.loadmat('slopes_ave_update20230305_bad.mat')  # slopes_ave bad
ngs_simu_slopes_ave_bad = ngs_simu_slopes_ave_bad['slopes_ave']
ngs_simu_xslope_ave_bad = ngs_simu_slopes_ave_bad[:, 0]  # bad
ngs_simu_yslope_ave_bad = ngs_simu_slopes_ave_bad[:, 1]  # bad
ngs_simu_slopes_ave1wfs_bad = scipy.io.loadmat('slopes_ave1wfs_update20230305_bad.mat')  # slopes_ave bad
ngs_simu_slopes_ave1wfs_bad = ngs_simu_slopes_ave1wfs_bad['slopes_ave1wfs']
ngs_simu_xslope_ave1wfs_bad = ngs_simu_slopes_ave1wfs_bad[:, 0]  # bad
ngs_simu_yslope_ave1wfs_bad = ngs_simu_slopes_ave1wfs_bad[:, 1]  # bad
ngs_simu_wave_tt_bad = scipy.io.loadmat('wave_tt_update20230305_bad.mat')  # slopes_ave bad
ngs_simu_wave_tt_bad = ngs_simu_wave_tt_bad['wave_tt']
ngs_simu_wave_tip_bad = ngs_simu_wave_tt_bad[:, 0]  # bad
ngs_simu_wave_tilt_bad = ngs_simu_wave_tt_bad[:, 1]  # bad



## 観測データの切り出し
ngs_time = ngs_obs[:, 0]*1e-3
ngs_simu_time = np.arange(1/400, (len(ngs_simu_xslope_ave_bad)+1)/400, 1/400)  # range(400Hz*time)

## 観測データの処理
ngs_obs_xslope_ave = ngs_obs[:, 1]  # obs [pix]
ngs_obs_yslope_ave = ngs_obs[:, 2]  # obs [pix]
ngs_obs_xslope_ave = ngs_obs_xslope_ave - np.mean(ngs_obs_xslope_ave)  # obs [pix]
ngs_obs_yslope_ave = ngs_obs_yslope_ave - np.mean(ngs_obs_yslope_ave)  # obs [pix]

## シミュレーションデータの処理
ngs_simu_xslope_ave_bad = ngs_simu_xslope_ave_bad - np.mean(ngs_simu_xslope_ave_bad)  # bad [pix]
ngs_simu_yslope_ave_bad = ngs_simu_yslope_ave_bad - np.mean(ngs_simu_yslope_ave_bad)  # bad [pix]
ngs_simu_xslope_ave1wfs_bad = ngs_simu_xslope_ave1wfs_bad - np.mean(ngs_simu_xslope_ave1wfs_bad)  # bad [arcsec]
ngs_simu_yslope_ave1wfs_bad = ngs_simu_yslope_ave1wfs_bad - np.mean(ngs_simu_yslope_ave1wfs_bad)  # bad [arcsec]
ngs_simu_wave_tip_bad = ngs_simu_wave_tip_bad - np.mean(ngs_simu_wave_tip_bad)  # bad [arcsec]
ngs_simu_wave_tilt_bad = ngs_simu_wave_tilt_bad - np.mean(ngs_simu_wave_tilt_bad)  # bad [arcsec]


## unit: [pix] -> [arcsec]
ngs_obs_xslope_ave = ngs_obs_xslope_ave*slope_unit_conversion  # obs [arcsec]
ngs_obs_yslope_ave = ngs_obs_yslope_ave*slope_unit_conversion  # obs [arcsec]

ngs_simu_xslope_ave_bad = ngs_simu_xslope_ave_bad*simu_unit_conversion  # bad [arcsec]
ngs_simu_yslope_ave_bad = ngs_simu_yslope_ave_bad*simu_unit_conversion  # bad [arcsec]


## パワースペクトル密度の計算
fs = 400  # サンプリング周波数[Hz]

pw_ngs_obs_xslope = signal.welch(ngs_obs_xslope_ave, fs)  # obs
pw_ngs_obs_yslope = signal.welch(ngs_obs_yslope_ave, fs)  # obs

pw_ngs_simu_xslope_bad = signal.welch(ngs_simu_xslope_ave_bad, fs)  # bad
pw_ngs_simu_yslope_bad = signal.welch(ngs_simu_yslope_ave_bad, fs)  # bad
pw_ngs_simu_xslope1wfs_bad = signal.welch(ngs_simu_xslope_ave1wfs_bad, fs)  # bad
pw_ngs_simu_yslope1wfs_bad = signal.welch(ngs_simu_yslope_ave1wfs_bad, fs)  # bad
pw_ngs_simu_wave_tip_bad = signal.welch(ngs_simu_wave_tip_bad, fs)  # bad
pw_ngs_simu_wave_tilt_bad = signal.welch(ngs_simu_wave_tilt_bad, fs)  # bad


## データのプロット
fig_1 = plt.figure(figsize=(9, 9))
ax1 = plt.subplot(3, 1, 1)
ax1.plot(ngs_time, ngs_obs_xslope_ave, color="k", linewidth=0.2, label="NGS_obs_xslope")
ax1.plot(ngs_simu_time, ngs_simu_xslope_ave_bad, color="green", linewidth=0.2, label="NGS_simu_xslope_bad")
ax1.plot(ngs_simu_time, ngs_simu_xslope_ave1wfs_bad, color="red", linewidth=0.2, label="NGS_simu_xslope1wfs_bad")
ax1.plot(ngs_simu_time, ngs_simu_wave_tip_bad, color="blue", linewidth=0.2, label="NGS_simu_wave_tip_bad")
ax1.set_ylabel("xslope [arcsec]")
ax1.set_xlabel("Time [s]")
ax1.set_xlim(0, 10)
ax1.legend()

ax2 = plt.subplot(3, 1, 2)
ax2.plot(ngs_time, ngs_obs_yslope_ave, color="k", linewidth=0.2, label="NGS_obs_yslope")
ax2.plot(ngs_simu_time, ngs_simu_yslope_ave_bad, color="green", linewidth=0.2, label="NGS_simu_yslope_bad")
ax2.plot(ngs_simu_time, ngs_simu_yslope_ave1wfs_bad, color="red", linewidth=0.2, label="NGS_simu_yslope1wfs_bad")
ax2.plot(ngs_simu_time, ngs_simu_wave_tilt_bad, color="blue", linewidth=0.2, label="NGS_simu_wave_tilt_bad")
ax2.set_ylabel("yslope [arcsec]")
ax2.set_xlabel("Time [s]")
ax2.set_xlim(0, 10)
ax2.legend()

# パワースペクトル密度のプロット
ax3 = plt.subplot(3, 1, 3)
ax3.plot(pw_ngs_obs_xslope[0], pw_ngs_obs_xslope[1], color="k", linewidth=1, label="NGS_obs_xslope")
ax3.plot(pw_ngs_obs_yslope[0], pw_ngs_obs_yslope[1], color="k", linewidth=1, linestyle="dashed", label="NGS_obs_yslope")
ax3.plot(pw_ngs_simu_xslope_bad[0], pw_ngs_simu_xslope_bad[1], color="green", linewidth=1, label="NGS_simu_xslope_bad")
ax3.plot(pw_ngs_simu_yslope_bad[0], pw_ngs_simu_yslope_bad[1], color="green", linewidth=1, linestyle="dashed", label="NGS_simu_yslope_bad")
ax3.plot(pw_ngs_simu_xslope1wfs_bad[0], pw_ngs_simu_xslope1wfs_bad[1], color="red", linewidth=1, label="NGS_simu_xslope1wfs_bad")
ax3.plot(pw_ngs_simu_yslope1wfs_bad[0], pw_ngs_simu_yslope1wfs_bad[1], color="red", linewidth=1, linestyle="dashed", label="NGS_simu_yslope1wfs_bad")
ax3.plot(pw_ngs_simu_wave_tip_bad[0], pw_ngs_simu_wave_tip_bad[1], color="blue", linewidth=1, label="NGS_simu_wave_tip_bad")
ax3.plot(pw_ngs_simu_wave_tilt_bad[0], pw_ngs_simu_wave_tilt_bad[1], color="blue", linewidth=1, linestyle="dashed", label="NGS_simu_wave_tilt_bad")
ax3.set_xlabel("Frequency [Hz]")
ax3.set_ylabel("PSD")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlim(1, 200)
#ax3.set_ylim(1e-6, 1e+2)
ax3.legend()

plt.show()