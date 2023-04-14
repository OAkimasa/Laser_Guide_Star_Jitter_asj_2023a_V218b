import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

from VectorFunctions import VectorFunctions

start = time.time()


N_air = 1.0
N_Fused_Silica = 1.458  # Refractive index (Fused Silica at 589nm)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect((1, 1, 1))
ax.set_xlim(0, 150)
ax.set_ylim(-75, 75)
ax.set_zlim(-75, 75)


## ---- 単位は mm ----

## param = [pos, nV, R, Lens_R(axis+-)]
surface_1 = [[0.000, 0., 0.], [1., 0., 0.], 25.4/2, 51.09]  # air->lens1
surface_2 = [[3.000, 0., 0.], [1., 0., 0.], 25.4/2, -46.70]  # lens1->air
surface_3 = [[98.693, 0., 0.], [1., 0., 0.], 38.1/2, 130.2]  # air->lens2
surface_4 = [[105.693, 0., 0.], [1., 0., 0.], 38.1/2, 67.95]  # lens2->air
surface_5 = [[148.000, 0., 0.], [1., 0., 0.], 50.8/2, 842.33]  # air->lens3
surface_6 = [[155.000, 0., 0.], [1., 0., 0.], 50.8/2, 145.34]  # lens3->air
evaluation_plane = [[200.000, 0., 0.], [1., 0., 0.], 0., 0.]  # air->air
surface_list = [surface_1, surface_2, surface_3, surface_4, surface_5, surface_6, evaluation_plane]


## インスタンス生成, 光線追跡開始
VF = VectorFunctions()
VF.ax = ax  # axをVFに登録


for surface in surface_list:
    VF.plot_lens(surface)  # レンズ描画


# 始点を生成する
width = 3
space = 1
rayDensity = 1
rayCenterX = -30
rayCenterY = 0
rayCenterZ = 0
size = len(np.arange(-width+rayCenterY, 1+width+rayCenterY, space))**2
pointsY, pointsZ = np.meshgrid(
    np.arange(-width+rayCenterY, 1+width+rayCenterY, space),
    np.arange(-width+rayCenterY, 1+width+rayCenterY, space))
pointsX = np.array([rayCenterX]*size)
pointsY = pointsY.reshape(size)*rayDensity
pointsZ = pointsZ.reshape(size)*rayDensity
raySPoint0 = VF.makePoints(pointsX, pointsY, pointsZ, size, 3)


# 初期値
#ray_start_pos_init = np.array([-30.0, 3.0, 3.0])  # 初期値
ray_start_pos_init = raySPoint0  # 初期値
ray_start_dir_init = np.array([[1.0, 0.0, 0.0]]*len(raySPoint0))  # 初期値

# surface_1
VF.ray_start_pos = ray_start_pos_init  # 初期値
VF.ray_start_dir = ray_start_dir_init  # 初期値
VF.set_surface(surface_1)  # surface_1をVFに登録
VF.raytrace_sphere()  # 光線追跡
VF.refractive_index_before = N_air  # 入射側の屈折率
VF.refractive_index_after = N_Fused_Silica  # 出射側の屈折率
VF.refract()  # 空気からレンズ1の屈折
VF.plotLineRed()  # 光線描画

# surface_2
VF.ray_start_pos = VF.ray_end_pos  # surface_1の終点をsurface_2の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_1の終点をsurface_2の始点に
VF.set_surface(surface_2)  # surface_2をVFに登録
VF.raytrace_sphere()  # 光線追跡
VF.refractive_index_before = N_Fused_Silica  # 入射側の屈折率
VF.refractive_index_after = N_air  # 出射側の屈折率
VF.refract()  # レンズ1から空気の屈折
VF.plotLineRed()  # 光線描画

# surface_3
VF.ray_start_pos = VF.ray_end_pos  # surface_2の終点をsurface_3の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_2の終点をsurface_3の始点に
VF.set_surface(surface_3)  # surface_3をVFに登録
VF.raytrace_sphere()  # 光線追跡
VF.refractive_index_before = N_air  # 入射側の屈折率
VF.refractive_index_after = N_Fused_Silica  # 出射側の屈折率
VF.refract()  # 空気からレンズ2の屈折
VF.plotLineRed()  # 光線描画

# surface_4
VF.ray_start_pos = VF.ray_end_pos  # surface_3の終点をsurface_4の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_3の終点をsurface_4の始点に
VF.set_surface(surface_4)  # surface_4をVFに登録
VF.raytrace_sphere()  # 光線追跡
VF.refractive_index_before = N_Fused_Silica  # 入射側の屈折率
VF.refractive_index_after = N_air  # 出射側の屈折率
VF.refract()  # レンズ2から空気の屈折
VF.plotLineRed()  # 光線描画

# surface_5
VF.ray_start_pos = VF.ray_end_pos  # surface_4の終点をsurface_5の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_4の終点をsurface_5の始点に
VF.set_surface(surface_5)  # surface_5をVFに登録
VF.raytrace_sphere()  # 光線追跡
VF.refractive_index_before = N_air  # 入射側の屈折率
VF.refractive_index_after = N_Fused_Silica  # 出射側の屈折率
VF.refract()  # 空気からレンズ3の屈折
VF.plotLineRed()  # 光線描画

# surface_6
VF.ray_start_pos = VF.ray_end_pos  # surface_5の終点をsurface_6の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_5の終点をsurface_6の始点に
VF.set_surface(surface_6)  # surface_6をVFに登録
VF.raytrace_sphere()  # 光線追跡
VF.refractive_index_before = N_Fused_Silica  # 入射側の屈折率
VF.refractive_index_after = N_air  # 出射側の屈折率
VF.refract()  # レンズ3から空気の屈折
VF.plotLineRed()  # 光線描画

# evaluation plane
VF.ray_start_pos = VF.ray_end_pos  # surface_6の終点を評価面の始点に
VF.ray_start_dir = VF.ray_end_dir  # surface_6の終点を評価面の始点に
VF.set_surface(evaluation_plane)  # 評価面をVFに登録
VF.raytrace_plane()  # 光線追跡
VF.refractive_index_before = N_air  # 入射側の屈折率
VF.refractive_index_after = N_air  # 出射側の屈折率
VF.plotLineRed()  # 光線描画

"""focal_length = VF.calcFocalLength(ray_start_pos_init)  # 計算した焦点距離を取得
print("焦点距離 = " + str(focal_length) + " mm")"""


print("run time: {0:.3f} sec".format(time.time() - start))
plt.show()