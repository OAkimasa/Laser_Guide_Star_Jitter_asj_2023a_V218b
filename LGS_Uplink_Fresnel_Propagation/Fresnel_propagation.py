import numpy as np
import math
import scipy.optimize as opt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.io
import aotools
import poppy
import astropy.units as u
import csv
import tqdm
import time


N_pix = 128 # number of pixels
nx_size = 1.0 # meter 
pixel_scale = nx_size / float(N_pix) # meter/pix
wl_laser = 589.0 # laser wavelength in nm
t_alt = np.array([0., 500., 1000., 2000., 4000., 8000., 16000.]) # Turbulence height in m
EL = np.array([30.0, 55.0, 60.0, 90.0]) # telescope elevation in deg
X = (1.0/np.cos(np.radians(90.0-EL))) # Airmass
unit_conversion = np.rad2deg(np.arcsin(1/(90*1000/np.cos(np.radians(90.0-55)))))*3600/128  # arcsec/pix

x,y = np.meshgrid(np.arange(0,N_pix,1), np.arange(0,N_pix,1))
R = np.sqrt((x-N_pix/2.)**2 + (y-N_pix/2.)**2)

w0 = 0.112 # m
w0_pix = w0 * math.sqrt(2.0) / pixel_scale # multiplyt by sqrt(2) to adjust the waist size in the fresnel propagator
gauss =np.exp(-2*R**2/w0_pix**2)

def make_LLT_mask(target_diameter_arcsec=40):
    # mask
    mask_x = np.arange(0, 128, 1)
    mask_y = np.arange(0, 128, 1)
    mask_X, mask_Y = np.meshgrid(mask_x, mask_y)

    def param_func(target_diameter_arcsec):
        """if target_diameter_arcsec == 10:
            r_2_x = 0.0685
            r_2_y = 0.0666
            mask4_line = 0.138
        elif target_diameter_arcsec == 20:
            r_2_x = 0.0628
            r_2_y = 0.0642
            mask4_line = 0.2043
        elif target_diameter_arcsec == 30:
            r_2_x = -0.0874
            r_2_y = -0.0882
            mask4_line = 0.1386
        elif target_diameter_arcsec == 40:
            r_2_x = -0.1025
            r_2_y = -0.1034
            mask4_line = 0.2  # tmp"""

        # r_2_xy (ビーム中心とLLT開口中心の距離[m]), 3次関数でフィッティング
        fitted_curve = (np.poly1d(np.polyfit([10, 20, 30, 40],
                                             [0.0685, 0.0628, -0.0874, -0.1025], 3)))
        r_2_x = fitted_curve(target_diameter_arcsec)
        r_2_y = fitted_curve(target_diameter_arcsec)

        # プリズムのケラレ位置とLLT開口中心の距離[m], 1次関数でフィッティング
        fitted_line_tmp = (np.poly1d(np.polyfit([9, 10, 15, 18, 20, 23, 28, 30, 32, 34], [
            -431, -445, -474, -490, -505, -526, -565, -596, -594, -615], 1))+369.6)*1e-3
        """fitted_line = np.poly1d(np.polyfit(
            [10, 20, 30], [0.138, 0.2043, 0.1386], 1))"""
        # プリズムのケラレ位置 - ビーム中心位置
        param_mask4_line = -1. * \
            (fitted_line_tmp(target_diameter_arcsec)-r_2_x)
        #print("r_2_x: ", r_2_x)
        #print("param_mask4_line: ", param_mask4_line)

        return [r_2_x, r_2_y, param_mask4_line]
    r_2_x, r_2_y, param_mask4_line = param_func(target_diameter_arcsec)
    mask4_line = N_pix/2 - param_mask4_line/(1./128.)

    # mask1
    # beam waistでカットするマスク
    r_1 = np.sqrt((mask_X-N_pix/2)**2 + (mask_Y-N_pix/2)**2)
    mask1_radius = 0.247  # m
    mask1 = np.copy(r_1)*0.0
    idx1_1 = np.where(r_1 >= mask1_radius / (1./128.))
    idx1_2 = np.where(r_1 < mask1_radius / (1./128.))
    mask1[idx1_1] = 0.0
    mask1[idx1_2] = 1.0

    # mask2
    # LLT開口マスク
    mask2_radius = 0.247  # m
    r_2 = np.sqrt((mask_X-N_pix/2+r_2_x*N_pix)**2 +
                  (mask_Y-N_pix/2+r_2_y*N_pix)**2)
    mask2 = np.copy(r_2)*0.0
    idx2_1 = np.where(r_2 >= mask2_radius / (1./128.))
    idx2_2 = np.where(r_2 < mask2_radius / (1./128.))
    mask2[idx2_1] = 0.0
    mask2[idx2_2] = 1.0

    # mask3
    # LLT副鏡マスク
    r_3 = r_2
    mask3_radius = 0.026  # m
    mask3 = np.copy(r_3)*0.0
    idx3_1 = np.where(r_3 >= mask3_radius / (1./128.))
    idx3_2 = np.where(r_3 < mask3_radius / (1./128.))
    mask3[idx3_1] = 1.0
    mask3[idx3_2] = 0.0

    # mask4
    # プリズムマスク
    mask4 = np.copy(mask_X)*0.0
    idx4_1 = np.where(mask_X >= mask4_line)
    idx4_2 = np.where(mask_X < mask4_line)
    mask4[idx4_1] = 1.0
    mask4[idx4_2] = 0.0

    # mask5
    # プリズムマスク
    mask5_line = mask4_line
    mask5 = np.copy(mask_Y)*0.0
    idx5_1 = np.where(mask_Y >= mask5_line)
    idx5_2 = np.where(mask_Y < mask5_line)
    mask5[idx5_1] = 1.0
    mask5[idx5_2] = 0.0

    image = gauss*mask1*mask2*mask3*mask4*mask5
    return image

def gauss2d(pos, A, x0, y0, w, q, pa):
    x = pos[0]
    y = pos[1]
    xc = (x - x0) * np.cos(np.deg2rad(-pa)) - (y - y0) * np.sin(np.deg2rad(-pa))
    yc = (x - x0) * np.sin(np.deg2rad(-pa)) + (y - y0) * np.cos(np.deg2rad(-pa))
    r = np.sqrt(xc**2 + (yc/q)**2)
    g = A * np.exp(-(2.0*r**2/(w**2)))
    return g.ravel()

# phase_screenはOOMAOで計算した。simulate_turbulence.mを参照
def load_phase_screen_200Hz_EL55_median():
    print("loading phase screens...")
    loadtime = time.time()
    # load phase screen data
    """dict_1 = scipy.io.loadmat('phase_200Hz_EL55_part1.mat')
    dict_2 = scipy.io.loadmat('phase_200Hz_EL55_part2.mat')
    dict_3 = scipy.io.loadmat('phase_200Hz_EL55_part3.mat')"""
    dict_1 = scipy.io.loadmat('phase_200Hz_EL55_UPDATE20230303_part1.mat')  # median
    dict_2 = scipy.io.loadmat('phase_200Hz_EL55_UPDATE20230303_part2.mat')  # median
    dict_3 = scipy.io.loadmat('phase_200Hz_EL55_UPDATE20230303_part3.mat')  # median

    phase_screens = np.concatenate((dict_1['phase'], dict_2['phase'], dict_3['phase']), axis=2)
    #print("phase_screens.shape = ", phase_screens.shape)

    # 6000, 7, 128, 128に変形
    phscrn_el = []
    for i in range(6000):
        phscrn_arr = []
        for j in range(7):
            phscrn_arr.append(phase_screens[:, :, i, j])
        phscrn_el.append(np.array(phscrn_arr))
    print("complete!")
    print("loadtime = ", time.time()-loadtime)
    return phscrn_el

def load_phase_screen_200Hz_EL55_good():
    print("loading phase screens...")
    loadtime = time.time()
    # load phase screen data
    dict_1 = scipy.io.loadmat('phase_200Hz_EL55_update20230305_good_part1.mat')  # good

    phase_screens = dict_1['phase']  # 1のみ
    #print("phase_screens.shape = ", phase_screens.shape)

    # 2000, 7, 128, 128に変形
    phscrn_el = []
    for i in range(2000):
        phscrn_arr = []
        for j in range(7):
            phscrn_arr.append(phase_screens[:, :, i, j])
        phscrn_el.append(np.array(phscrn_arr))
    print("complete!")
    print("loadtime = ", time.time()-loadtime)
    return phscrn_el

def load_phase_screen_200Hz_EL55_bad():
    print("loading phase screens...")
    loadtime = time.time()
    # load phase screen data
    dict_1 = scipy.io.loadmat('phase_200Hz_EL55_update20230305_bad_part1.mat')  # bad

    phase_screens = dict_1['phase']  # 1のみ
    #print("phase_screens.shape = ", phase_screens.shape)

    # 2000, 7, 128, 128に変形
    phscrn_el = []
    for i in range(2000):
        phscrn_arr = []
        for j in range(7):
            phscrn_arr.append(phase_screens[:, :, i, j])
        phscrn_el.append(np.array(phscrn_arr))
    print("complete!")
    print("loadtime = ", time.time()-loadtime)
    return phscrn_el

#start_image = gauss  # no mask
#start_image = make_beam_image_40arcsec()  # LLT mask
target_diameter_arcsec = 18
start_image = make_LLT_mask(target_diameter_arcsec=target_diameter_arcsec)  # LLT mask
print("start_image_sum = ", np.sum(start_image))
"""fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
Z = start_image.reshape(N_pix, N_pix)
ax.plot_surface(x, y, Z)"""
#plt.show()


roop_num = 1     # 1 frame
#roop_num = 2000  # frames
popt_list = []
total_photon_list = []
fit_total_photon_list = []
res_total_photon_list = []
LGS_intensity_list = []
flag_phase_screen = 1 # 0: No phase screen, 1: phase screen


#phscrn_el_list = load_phase_screen_200Hz_EL55_median()  # OOMAO,ファイルから読み込み
#phscrn_el_list = load_phase_screen_200Hz_EL55_good()  # OOMAO,ファイルから読み込み
phscrn_el_list = load_phase_screen_200Hz_EL55_bad()  # OOMAO,ファイルから読み込み
#print("phscrn_el_list.shape = ", np.array(phscrn_el_list).shape)

for i in tqdm.tqdm(range(roop_num)):
    i_EL = 1  # Airmass index  EL30:0, EL55:1, EL60:2, EL90:3 -> Airmass:2.0, 1.221, 1.1547..., 1.0
    #i = 5999                           # select 1 frame
    N_frame = str(i+1).zfill(5)
    print('\nPhase screen number = ', N_frame)
    phscrn_el = phscrn_el_list[i]  # select time
    #print("phscrn_el.shape = ", np.array(phscrn_el).shape)


    # Setup initial condition for the fresnel propagation
    wf = poppy.FresnelWavefront(0.5*u.m, wavelength=wl_laser*1.0e-9, npix=N_pix, oversample=4)
    transmission = start_image*0.0+1.0
    #wf.display()
    #plt.show()

    opd = phscrn_el[0]*(wl_laser*1.0e-9/(2.*math.pi))*flag_phase_screen
    input_wave = poppy.ArrayOpticalElement(transmission=start_image, opd=opd*u.m, name='input wavefront', pixelscale=pixel_scale*u.m/u.pixel)
    wf *= input_wave
    #wf.display(what='both')
    #plt.show()

    # propagation
    wf.propagate_fresnel(0.5e+3*X[i_EL]*u.m)
    opd = phscrn_el[1]*(wl_laser*1.0e-9/(2.*math.pi))*flag_phase_screen
    ph_05 = poppy.ArrayOpticalElement(transmission=transmission, opd=opd*u.m, name='0.5km', pixelscale=pixel_scale*u.m/u.pixel)
    wf *= ph_05
    #wf.display(what='both')
    #plt.show()

    wf.propagate_fresnel(0.5e+3*X[i_EL]*u.m)
    opd = phscrn_el[2]*(wl_laser*1.0e-9/(2.*math.pi))*flag_phase_screen
    ph_1 = poppy.ArrayOpticalElement(transmission=transmission, opd=opd*u.m, name='1km', pixelscale=pixel_scale*u.m/u.pixel)
    wf *= ph_1
    #wf.display(what='both')
    #plt.show()

    wf.propagate_fresnel(1.e+3*X[i_EL]*u.m)
    opd = phscrn_el[3]*(wl_laser*1.0e-9/(2.*math.pi))*flag_phase_screen
    ph_2 = poppy.ArrayOpticalElement(transmission=transmission, opd=opd*u.m, name='2km', pixelscale=pixel_scale*u.m/u.pixel)
    wf *= ph_2
    #wf.display(what='both')
    #plt.show()

    wf.propagate_fresnel(2.e+3*X[i_EL]*u.m)
    opd = phscrn_el[4]*(wl_laser*1.0e-9/(2.*math.pi))*flag_phase_screen
    ph_4 = poppy.ArrayOpticalElement(transmission=transmission, opd=opd*u.m, name='4km', pixelscale=pixel_scale*u.m/u.pixel)
    wf *= ph_4
    #wf.display(what='both')
    #plt.show()

    wf.propagate_fresnel(4.e+3*X[i_EL]*u.m)
    opd = phscrn_el[5]*(wl_laser*1.0e-9/(2.*math.pi))*flag_phase_screen
    ph_8 = poppy.ArrayOpticalElement(transmission=transmission, opd=opd*u.m, name='8km', pixelscale=pixel_scale*u.m/u.pixel)
    wf *= ph_8
    #wf.display(what='both')
    #plt.show()

    wf.propagate_fresnel(8.e+3*X[i_EL]*u.m)
    opd = phscrn_el[6]*(wl_laser*1.0e-9/(2.*math.pi))*flag_phase_screen
    ph_16 = poppy.ArrayOpticalElement(transmission=transmission, opd=opd*u.m, name='16km', pixelscale=pixel_scale*u.m/u.pixel)
    wf *= ph_16
    #wf.display(what='both')
    #plt.show()

    # last propagation
    wf.propagate_fresnel((90.e+3-16.e+3)*X[i_EL]*u.m)
    #wf.display(what='both')
    N_frame = str(i+1).zfill(5)
    #plt.savefig('frame_'+ N_frame +'.png')
    #plt.show()


    # 2D gaussian fit
    # A, x0, y0, w, q, pa
    gauss_peak = np.max(wf.intensity)
    initial_guess = (gauss_peak, 2.0*N_pix, 2.0*N_pix, 0.1*w0/pixel_scale, 1.0, 0)
    Xp,Yp = np.meshgrid(np.arange(0,N_pix*4,1), np.arange(0,N_pix*4,1))
    popt, pcov = opt.curve_fit(gauss2d, (Xp, Yp), wf.intensity.ravel(), p0=initial_guess)
    print('popt :', popt)
    popt_list.append(popt)


    Z_initial_guess = gauss2d((Xp,Yp), *initial_guess).reshape(Xp.shape)
    Z_intensity = wf.intensity.reshape(Xp.shape)
    LGS_intensity_list.append(Z_intensity)
    #Z_fit = gauss2d((Xp,Yp), *popt).reshape(Xp.shape)

    """total_photon = np.sum(Z_intensity)
    print("total_photon : ", total_photon)
    total_photon_list.append([total_photon])

    fit_total_photon = np.sum(Z_fit)
    print("fit_photon : ", fit_total_photon)
    fit_total_photon_list.append([fit_total_photon])

    res_total_photon = np.sum(Z_intensity-Z_fit)
    print("res_total_photon : ", res_total_photon)
    res_total_photon_list.append([res_total_photon])"""

    fig = plt.figure(figsize=(9,9))
    plt.subplot(1,1,1)
    plt.imshow(Z_intensity, cmap='jet', vmax=0.35, vmin=0.0)
    plt.colorbar()
    #plt.savefig('frame_'+ N_frame +'.png')
    plt.show()

    """fig = plt.figure(figsize=(9,9))
    plt.subplot(1,1,1)
    plt.imshow(Z_fit, cmap='jet', vmax=0.35, vmin=0.0)
    plt.colorbar()
    #plt.savefig('frame_'+ N_frame +'.png')
    plt.show()"""

    """fig = plt.figure(figsize=(9,9))
    plt.subplot(1,1,1)
    plt.imshow(Z_intensity-Z_fit, cmap='jet', vmax=0.35, vmin=0.0)
    plt.colorbar()
    plt.savefig('frame_'+ N_frame +'.png')
    #plt.show()"""

    """fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(Xp, Yp, Z_intensity)
    ax.plot_wireframe(Xp, Yp, Z_fit, color='r')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_zlim(0, 0.3)
    ax.view_init(0, 45)
    #plt.savefig('frame_'+ N_frame +'.png')
    plt.show()"""

# save mat file
#scipy.io.savemat('LGS_intensity_200Hz_EL55_update_bad.mat', mdict={'intensity': LGS_intensity_list})

#print('popt_list :', popt_list)
"""with open('popt_list_200Hz_EL55_update_bad.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(popt_list)"""

#print('total_photon_list :', total_photon_list)
"""with open('total_photon_list.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(total_photon_list)"""

#print('fit_total_photon_list :', fit_total_photon_list)
"""with open('fit_total_photon_list.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(fit_total_photon_list)"""

#print('res_total_photon_list :', res_total_photon_list)
"""with open('res_total_photon_list.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(res_total_photon_list)"""