import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

from VectorFunctions import VectorFunctions


N_air = 1.0

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect((1, 1, 1))
ax.set_xlim(0, 50)
ax.set_ylim(-25, 25)
ax.set_zlim(-25, 25)


## ---- 単位は mm ----

## param = [pos, nV, R, Lens_R(axis+-)]
surface_1 = [[0., 0., 0.], [1., 0., 0.], 17.91, -47.8]  # air->lens1
surface_2 = [[7.56, 0., 0.], [1., 0., 0.], 17.91, 93.0]  # lens1->lens2
surface_3 = [[10.26, 0., 0.], [1., 0., 0.], 17.88, -262.5]  # lens2->air
surface_4 = [[18.86, 0., 0.], [1., 0., 0.], 12.60, 112.7]  # air->lens3
surface_5 = [[24.0, 0., 0.], [1., 0., 0.], 12.60, -43.0]  # lens3->air
surface_6 = [[36.76, 0., 0.], [1., 0., 0.], 13.95, 10000000]  # air->lens4
surface_7 = [[38.76, 0., 0.], [1., 0., 0.], 13.95, -65.8]  # lens4->lens5
surface_8 = [[44.56, 0., 0.], [1., 0., 0.], 13.95, 93.05]  # lens5->air
surface_list = [surface_1, surface_2, surface_3, surface_4, surface_5, surface_6, surface_7, surface_8]

# glass_list: [Fraunhofer C, Fraunhofer F]
N_lens_1 = [1.60008, 1.61003]  # N-SK14
N_lens_2 = [1.49552, 1.50296]  # N-BK10
N_lens_3 = [1.61506, 1.63208]  # N-F2
N_lens_4 = [1.50592, 1.51423]  # N-ZK7
N_lens_5 = [1.58619, 1.59581]   # N-SK5
paramList = [N_lens_1, N_lens_2, N_lens_3, N_lens_4, N_lens_5]


## レンズ描画
def plot_lens(params):
    geneNum = 300
    limitTheta = 2*np.pi  # theta生成数
    limitPhi = np.pi  # phi生成数
    theta = np.linspace(0, limitTheta, geneNum)
    phi = np.linspace(0, limitPhi, geneNum)

    argmax_index = np.argmax(np.abs(params[1]))

    if argmax_index == 0:
        Ys = np.outer(np.sin(theta), np.sin(phi))
        Zs = np.outer(np.ones(np.size(theta)), np.cos(phi))
        Ys1 = params[2] * Ys
        Zs1 = params[2] * Zs
        if params[3] < 0:
            Xs1 = -(params[3]**2-Ys1**2-Zs1**2)**0.5 - params[3]
            ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)
        elif params[3] > 0:
            Xs1 = (params[3]**2-Ys1**2-Zs1**2)**0.5 - params[3]
            ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)
    elif argmax_index == 1:
        Xs = np.outer(np.sin(theta), np.sin(phi))
        Zs = np.outer(np.ones(np.size(theta)), np.cos(phi))
        Xs1 = params[2] * Xs
        Zs1 = params[2] * Zs
        if params[3] < 0:
            Ys1 = -(params[3]**2-Xs1**2-Zs1**2)**0.5 - params[3]
            ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)
        elif params[3] > 0:
            Ys1 = (params[3]**2-Xs1**2-Zs1**2)**0.5 - params[3]
            ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)
    elif argmax_index == 2:
        Xs = np.outer(np.sin(theta), np.sin(phi))
        Ys = np.outer(np.ones(np.size(theta)), np.cos(phi))
        Xs1 = params[2] * Xs
        Ys1 = params[2] * Ys
        if params[3] < 0:
            Zs1 = -(params[3]**2-Xs1**2-Ys1**2)**0.5 - params[3]
            ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)
        elif params[3] > 0:
            Zs1 = (params[3]**2-Xs1**2-Ys1**2)**0.5 - params[3]
            ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)

for surface in surface_list:
    plot_lens(surface)  # レンズ描画


## インスタンス生成, 光線追跡開始
VF = VectorFunctions()
VF.ax = ax  # axをVFに登録

# lambda_index, 0:Fraunhofer C, 1:Fraunhofer F
lambda_index = 0  # Fraunhofer C
#lambda_index = 1  # Fraunhofer F

# 初期値
ray_start_pos_init = np.array([-5.0, 8.0, 8.0])  # 初期値
ray_start_dir_init = np.array([1.0, 0.0, 0.0])  # 初期値

# surface_1
VF.ray_start_pos = ray_start_pos_init  # 初期値
VF.ray_start_dir = ray_start_dir_init  # 初期値
pos_surface_1 = np.array(surface_1[0])  # [x, y, z]
lens_R_surface_1 = surface_1[3]  # lens_R
VF.raytrace_sphere(pos_surface_1, lens_R_surface_1)  # 光線追跡
nV_surface_1 = VF.calcNormalV_sphere(pos_surface_1, lens_R_surface_1)  # 法線ベクトル
VF.refract(nV_surface_1, N_air, N_lens_1[lambda_index])  # 空気からレンズ1の屈折
VF.plotLineRed()  # 光線描画

# surface_2
VF.ray_start_pos = VF.ray_end_pos  # surface_1の終点をsurface_2の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_1の終点をsurface_2の始点に
pos_surface_2 = np.array(surface_2[0])  # [x, y, z]
lens_R_surface_2 = surface_2[3]  # lens_R
VF.raytrace_sphere(pos_surface_2, lens_R_surface_2)  # 光線追跡
nV_surface_2 = VF.calcNormalV_sphere(pos_surface_2, lens_R_surface_2)  # 法線ベクトル
VF.refract(nV_surface_2, N_lens_1[lambda_index], N_lens_2[lambda_index])  # レンズ1からレンズ2の屈折
VF.plotLineRed()  # 光線描画

# surface_3
VF.ray_start_pos = VF.ray_end_pos  # surface_2の終点をsurface_3の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_2の終点をsurface_3の始点に
pos_surface_3 = np.array(surface_3[0])  # [x, y, z]
lens_R_surface_3 = surface_3[3]  # lens_R
VF.raytrace_sphere(pos_surface_3, lens_R_surface_3)  # 光線追跡
nV_surface_3 = VF.calcNormalV_sphere(pos_surface_3, lens_R_surface_3)  # 法線ベクトル
VF.refract(nV_surface_3, N_lens_2[lambda_index], N_air)  # レンズ2から空気の屈折
VF.plotLineRed()  # 光線描画

# surface_4
VF.ray_start_pos = VF.ray_end_pos  # surface_3の終点をsurface_4の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_3の終点をsurface_4の始点に
pos_surface_4 = np.array(surface_4[0])  # [x, y, z]
lens_R_surface_4 = surface_4[3]  # lens_R
VF.raytrace_sphere(pos_surface_4, lens_R_surface_4)  # 光線追跡
nV_surface_4 = VF.calcNormalV_sphere(pos_surface_4, lens_R_surface_4)  # 法線ベクトル
VF.refract(nV_surface_4, N_air, N_lens_3[lambda_index])  # 空気からレンズ3の屈折
VF.plotLineRed()  # 光線描画

# surface_5
VF.ray_start_pos = VF.ray_end_pos  # surface_4の終点をsurface_5の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_4の終点をsurface_5の始点に
pos_surface_5 = np.array(surface_5[0])  # [x, y, z]
lens_R_surface_5 = surface_5[3]  # lens_R
VF.raytrace_sphere(pos_surface_5, lens_R_surface_5)  # 光線追跡
nV_surface_5 = VF.calcNormalV_sphere(pos_surface_5, lens_R_surface_5)  # 法線ベクトル
VF.refract(nV_surface_5, N_lens_3[lambda_index], N_air)  # レンズ3から空気の屈折
VF.plotLineRed()  # 光線描画

# surface_6
VF.ray_start_pos = VF.ray_end_pos  # surface_5の終点をsurface_6の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_5の終点をsurface_6の始点に
pos_surface_6 = np.array(surface_6[0])  # [x, y, z]
lens_R_surface_6 = surface_6[3]  # lens_R
VF.raytrace_sphere(pos_surface_6, lens_R_surface_6)  # 光線追跡
nV_surface_6 = VF.calcNormalV_sphere(pos_surface_6, lens_R_surface_6)  # 法線ベクトル
VF.refract(nV_surface_6, N_air, N_lens_4[lambda_index])  # 空気からレンズ4の屈折
VF.plotLineRed()  # 光線描画

# surface_7
VF.ray_start_pos = VF.ray_end_pos  # surface_6の終点をsurface_7の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_6の終点をsurface_7の始点に
pos_surface_7 = np.array(surface_7[0])  # [x, y, z]
lens_R_surface_7 = surface_7[3]  # lens_R
VF.raytrace_sphere(pos_surface_7, lens_R_surface_7)  # 光線追跡
nV_surface_7 = VF.calcNormalV_sphere(pos_surface_7, lens_R_surface_7)  # 法線ベクトル
VF.refract(nV_surface_7, N_lens_4[lambda_index], N_lens_5[lambda_index])  # レンズ4からレンズ5の屈折
VF.plotLineRed()  # 光線描画

# surface_8
VF.ray_start_pos = VF.ray_end_pos  # surface_7の終点をsurface_8の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_7の終点をsurface_8の始点に
pos_surface_8 = np.array(surface_8[0])  # [x, y, z]
lens_R_surface_8 = surface_8[3]  # lens_R
VF.raytrace_sphere(pos_surface_8, lens_R_surface_8)  # 光線追跡
nV_surface_8 = VF.calcNormalV_sphere(pos_surface_8, lens_R_surface_8)  # 法線ベクトル
VF.refract(nV_surface_8, N_lens_5[lambda_index], N_air)  # レンズ5から空気の屈折
VF.plotLineRed()  # 光線描画

## 焦点距離の計算
"""focal_length = VF.calcFocalLength(ray_start_pos_init)  # 計算した焦点距離を取得
print("焦点距離 = " + str(focal_length) + " mm")"""

## 焦点位置の計算
focal_pos = VF.calcFocalPos(ray_start_pos_init)  # 計算した焦点位置を取得
print("焦点位置 = " + str(focal_pos) + " mm")

# evaluation plane
VF.ray_start_pos = VF.ray_end_pos  # surface_8の終点を評価面の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_8の終点を評価面の始点に
pos_plane = np.array([200, 0, 0])  # [x, y, z]
nV_eval_plane = np.array([1, 0, 0])  # 法線ベクトル
VF.raytrace_plane(pos_plane, nV_eval_plane)  # 光線追跡
VF.plotLineRed()  # 光線描画

focal_length = VF.calcFocalLength(ray_start_pos_init)  # 計算した焦点距離を取得
print("焦点距離 = " + str(focal_length) + " mm")


plt.show()