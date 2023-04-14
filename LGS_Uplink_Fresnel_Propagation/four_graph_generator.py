import numpy as np
import math
import scipy.optimize as opt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scipy.io
from scipy import interpolate
import aotools
import poppy
import astropy.units as u
import csv
import tqdm
import time


#unit_conversion = np.rad2deg(np.arcsin(1/(90*1000/np.cos(np.radians(90.0-55)))))*3600/128  # ~0.0147 arcsec/pix
unit_conversion = 1./128.  # 1/128 m/pix

def load_phase_screen_200Hz_EL55():
    print("loading phase screens...")
    loadtime = time.time()
    # load phase screen data
    dict1 = scipy.io.loadmat('phase_200Hz_EL55_update20230303_part1.mat')
    dict2 = scipy.io.loadmat('phase_200Hz_EL55_update20230303_part2.mat')
    dict3 = scipy.io.loadmat('phase_200Hz_EL55_update20230303_part3.mat')
    phase_screens = np.concatenate((dict1['phase'], dict2['phase'], dict3['phase']), axis=2)  # 結合
    #print("phase_screens.shape = ", phase_screens.shape)

    # 2000, 7, 128, 128に変形
    phscrn_el = []
    for i in range(6000):
        phscrn_arr = []
        for j in range(7):
            phscrn_arr.append(phase_screens[:, :, i, j])
        phscrn_el.append(np.array(phscrn_arr))
    print("complete!")
    print("loadtime = ", time.time()-loadtime)
    return phscrn_el

def gauss2d(pos, A, x0, y0, w, q, pa):
    x = pos[0]
    y = pos[1]
    xc = (x - x0) * np.cos(np.deg2rad(-pa)) - (y - y0) * np.sin(np.deg2rad(-pa))
    yc = (x - x0) * np.sin(np.deg2rad(-pa)) + (y - y0) * np.cos(np.deg2rad(-pa))
    r = np.sqrt(xc**2 + (yc/q)**2)
    g = A * np.exp(-(2.0*r**2/(w**2)))
    return g.ravel()


# phase screen, A-w graph, LGS intensity(2D), fitting mesh(3D)
roop_num = 1
#roop_num = 2000*3

# ax1
phscrn_el_list = load_phase_screen_200Hz_EL55()  # OOMAO,ファイルから読み込み
# ax3
dict1 = scipy.io.loadmat('LGS_intensity_200Hz_EL55_update_part1.mat')
dict2 = scipy.io.loadmat('LGS_intensity_200Hz_EL55_update_part2.mat')
dict3 = scipy.io.loadmat('LGS_intensity_200Hz_EL55_update_part3.mat')
LGS_intensity_list = np.concatenate((dict1['intensity'], dict2['intensity'], dict3['intensity']), axis=0)

# popt_list for A-w graph--------------------------------------------------------
popt_list = []
with open('popt_list_200Hz_EL55_update_1.csv', 'r') as f:
    reader = csv.reader(f)
    popt_list_1 = [row for row in reader]
popt_list_1 = np.array(popt_list_1)
popt_list_1 = popt_list_1.astype(np.float64)

with open('popt_list_200Hz_EL55_update_2.csv', 'r') as f:
    reader = csv.reader(f)
    popt_list_2 = [row for row in reader]
popt_list_2 = np.array(popt_list_2)
popt_list_2 = popt_list_2.astype(np.float64)

with open('popt_list_200Hz_EL55_update_3.csv', 'r') as f:
    reader = csv.reader(f)
    popt_list_3 = [row for row in reader]
popt_list_3 = np.array(popt_list_3)
popt_list_3 = popt_list_3.astype(np.float64)

# (2000, 6)*3 -> (6000, 6)
popt_list = np.concatenate([popt_list_1, popt_list_2, popt_list_3], axis=0)
popt_list = np.array(popt_list)

#A_list = []
I_list = []  # W/m^2 Holzlohner et al. 2010
w_list = []
print("calculating I_list...")
start_time = time.time()
for i in tqdm.tqdm(range(6000)):
    #A_list.append(popt_list[i][0])
    # IP2 ===============================
    P1 = 4.0
    pix_scale = 1.0/128.0 # m/pix
    total = np.sum(LGS_intensity_list[i])
    wf_n = LGS_intensity_list[i] / total
    I_peak = np.max(wf_n)
    IA = np.linspace(0,I_peak,100)
    fla = []
    for I in IA: 
        idx = np.where(wf_n > I)
        fla.append(np.sum(wf_n[idx]))
    ip = interpolate.interp1d(fla, IA)
    IP2 = ip(0.5)
    I = IP2*P1/pix_scale**2
    #print(I, 'W/m^2')
    I_list.append(I)
    # ===================================
    if abs(popt_list[i][4])<1:
        w_list.append(popt_list[i][3]*unit_conversion)
    else:
        w_list.append(abs(popt_list[i][4])*popt_list[i][3]*unit_conversion)
print("complete!")
print("time = ", time.time()-start_time)
#--------------------------------------------------------------------------


for i in tqdm.tqdm(range(roop_num)):
    #for i in tqdm.tqdm(np.array(A_upper_list)-1):
    #i = 1111                           # select 1 frame
    N_frame = str(i+1).zfill(5)
    #print('\nPhase screen number = ', N_frame)
    fig = plt.figure(figsize=(9,9))

    # Load the phase screen mat file
    ax1 = fig.add_subplot(2,2,1)
    phscrn_el = phscrn_el_list[i]  # select time
    phscrn_el = np.array(phscrn_el)  # list -> array, (7, 128, 128)
    phscrn_el = phscrn_el[0] + phscrn_el[1] + phscrn_el[2] + phscrn_el[3] + phscrn_el[4] + phscrn_el[5] + phscrn_el[6]  # (128, 128)
    plt.imshow(phscrn_el, cmap='jet')
    ax1.axis('off')
    cax = inset_axes(
        ax1,
        width="5%",
        height="90%",
        loc='lower left',
        bbox_to_anchor=(-0.1, 0.025, 1, 1),
        bbox_transform=ax1.transAxes,
        )  # left-align the colorbar
    plt.colorbar(cax=cax)
    cax.yaxis.set_ticks_position('left')
    ax1.set_title('phase screen '+N_frame, fontsize=9)

    # Load the popt_list_1~3 csv file
    ax2 = fig.add_subplot(2,2,2)
    plt.scatter(w_list, I_list, s=1, c='b', label='EL55')
    plt.scatter(0.212022, 24.922797, c='r', marker="x", label='without phase screen')
    plt.scatter(0.215337, 27.425417, c='r', marker="x", label='without phase screen & LLT mask')
    #plt.scatter(w_list[i], I_list[i], c='orange', marker="^", label='EL55, frame '+N_frame)
    #plt.xlabel('beam waist [m]', fontsize=10)
    #plt.ylabel('I [W/m^2]', fontsize=10)
    plt.xlim(0, 1.0)
    plt.ylim(0, 30)
    plt.grid()
    plt.legend(fontsize=7, loc='upper right')
    #ax2.set_title('I - beam waist graph', fontsize=9)

    # Load the LGS intensity mat file
    ax3 = fig.add_subplot(2,2,3)
    #print("LGS_intensity_list.shape = ", LGS_intensity_list.shape)
    Z_intensity = LGS_intensity_list[i]
    plt.imshow(Z_intensity, cmap='jet', vmax=0.30, vmin=0.0)
    ax3.axis('off')
    cax = inset_axes(
        ax3,
        width="5%",
        height="90%",
        loc='lower left',
        bbox_to_anchor=(-0.1, 0.025, 1, 1),
        bbox_transform=ax3.transAxes,
        )  # left-align the colorbar
    plt.colorbar(cax=cax)
    cax.yaxis.set_ticks_position('left')
    ax3.set_title('LGS at 90km EL55deg', fontsize=9)

    # Load the fitting mesh mat file
    ax4 = fig.add_subplot(2,2,4)
    fit_param = popt_list[i]
    #print("fit_param = ", fit_param)
    Xp, Yp = np.meshgrid(np.arange(0, 128*4, 1), np.arange(0, 128*4, 1))
    intensity_fit = gauss2d((Xp, Yp), *fit_param).reshape(Xp.shape)
    plt.imshow(intensity_fit, cmap='jet', vmax=0.30, vmin=0.0)
    ax4.axis('off')
    cax = inset_axes(
        ax4,
        width="5%",
        height="90%",
        loc='lower left',
        bbox_to_anchor=(-0.1, 0.025, 1, 1),
        bbox_transform=ax4.transAxes,
        )  # left-align the colorbar
    plt.colorbar(cax=cax)
    cax.yaxis.set_ticks_position('left')
    ax4.set_title('fitting', fontsize=9)

    #plt.savefig('frame_'+ N_frame +'_update20230306.png')
    plt.show()