import numpy as np

# 反射を計算する class
class VectorFunctions:
    def __init__(self):
        self.ax = None  # plotline関数のため
        self.ray_start_pos = np.array([0, 0, 0])  # 光線の始点
        self.ray_start_dir = np.array([0, 0, 0])  # 光線の方向ベクトル
        self.ray_end_pos = np.array([0, 0, 0])  # 光線の終点
        self.ray_end_dir = np.array([0, 0, 0])  # 光線の方向ベクトル
        self.surface_pos = np.array([0, 0, 0])  # 光学素子面中心の位置ベクトル
        self.lens_or_parabola_R = 0  # レンズの曲率半径またはパラボラの焦点距離
        self.refractive_index_before = 1  # 光線の入射側の屈折率
        self.refractive_index_after = 1  # 光線の出射側の屈折率
        self.normalV_optical_element = np.array([0, 0, 0])  # 光学素子中心の法線ベクトル
        self.normalV_refract_or_reflect = np.array([1, 0, 0])  # 屈折または反射計算に用いる法線ベクトル
        self.optical_path = 0  # 光路長


    ## 光学素子描画関数
    # ミラー描画
    def plot_mirror(self, params, color='gray'):
        X_center = params[0][0]
        Y_center = params[0][1]
        Z_center = params[0][2]
        normalV = params[1]
        R = params[2]

        argmax_index = np.argmax(np.abs(params[1]))

        if argmax_index == 0:
            # 円盤生成
            geneNum = 200
            y = np.linspace(Y_center-R, Y_center+R, geneNum)
            z = np.linspace(Z_center-R, Z_center+R, geneNum)
            Y, Z = np.meshgrid(y, z)
            if normalV[0] == 0:
                X = X_center - (normalV[1]*(Y-Y_center) +
                                normalV[2]*(Z-Z_center)) / 0.01
            else:
                X = X_center - (normalV[1]*(Y-Y_center) +
                                normalV[2]*(Z-Z_center)) / normalV[0]
            for i in range(geneNum):
                for j in range(geneNum):
                    if (X[i][j]-X_center)**2 + (Y[i][j]-Y_center)**2 + (Z[i][j]-Z_center)**2 > R**2:
                        Z[i][j] = np.nan

            #self.ax.quiver(X_center, Y_center, Z_center, normalV[0], normalV[1], normalV[2], color='black', length=50)
            self.ax.plot_wireframe(X, Y, Z, color=color, linewidth=0.3)
        elif argmax_index == 1:
            # 円盤生成
            geneNum = 200
            x = np.linspace(X_center-R, X_center+R, geneNum)
            z = np.linspace(Z_center-R, Z_center+R, geneNum)
            X, Z = np.meshgrid(x, z)
            if normalV[1] == 0:
                Y = Y_center - (normalV[0]*(X-X_center) +
                                normalV[2]*(Z-Z_center)) / 0.01
            else:
                Y = Y_center - (normalV[0]*(X-X_center) +
                                normalV[2]*(Z-Z_center)) / normalV[1]
            for i in range(geneNum):
                for j in range(geneNum):
                    if (X[i][j]-X_center)**2 + (Y[i][j]-Y_center)**2 + (Z[i][j]-Z_center)**2 > R**2:
                        Z[i][j] = np.nan

            #self.ax.quiver(X_center, Y_center, Z_center, normalV[0], normalV[1], normalV[2], color='black', length=50)
            self.ax.plot_wireframe(X, Y, Z, color=color, linewidth=0.3)
        elif argmax_index == 2:
            # 円盤生成
            geneNum = 200
            x = np.linspace(X_center-R, X_center+R, geneNum)
            y = np.linspace(Y_center-R, Y_center+R, geneNum)
            X, Y = np.meshgrid(x, y)
            if normalV[2] == 0:
                Z = Z_center - (normalV[0]*(X-X_center) +
                                normalV[1]*(Y-Y_center)) / 0.01
            else:
                Z = Z_center - (normalV[0]*(X-X_center) +
                                normalV[1]*(Y-Y_center)) / normalV[2]
            for i in range(geneNum):
                for j in range(geneNum):
                    if (X[i][j]-X_center)**2 + (Y[i][j]-Y_center)**2 + (Z[i][j]-Z_center)**2 > R**2:
                        Z[i][j] = np.nan

            #self.ax.quiver(X_center, Y_center, Z_center, normalV[0], normalV[1], normalV[2], color='black', length=50)
            self.ax.plot_wireframe(X, Y, Z, color=color, linewidth=0.3)

    def plot_window(self, params):
        X_center = params[0][0]
        Y_center = params[0][1]
        Z_center = params[0][2]
        normalV = params[1]
        R = params[2]

        argmax_index = np.argmax(np.abs(params[1]))

        if argmax_index == 0:
            # 円盤生成
            geneNum = 200
            y = np.linspace(Y_center-R, Y_center+R, geneNum)
            z = np.linspace(Z_center-R, Z_center+R, geneNum)
            Y, Z = np.meshgrid(y, z)
            if normalV[0] == 0:
                X = X_center - (normalV[1]*(Y-Y_center) +
                                normalV[2]*(Z-Z_center)) / 0.01
            else:
                X = X_center - (normalV[1]*(Y-Y_center) +
                                normalV[2]*(Z-Z_center)) / normalV[0]
            for i in range(geneNum):
                for j in range(geneNum):
                    if (X[i][j]-X_center)**2 + (Y[i][j]-Y_center)**2 + (Z[i][j]-Z_center)**2 > R**2:
                        Z[i][j] = np.nan

            #self.ax.quiver(X_center, Y_center, Z_center, normalV[0], normalV[1], normalV[2], color='black', length=50)
            self.ax.plot_wireframe(X, Y, Z, color='lightcyan', linewidth=0.3)
            theta = np.linspace(0, 2*np.pi, 100)
            x = np.zeros_like(theta) + X_center
            y = np.cos(theta)*R + Y_center
            z = np.sin(theta)*R + Z_center
            self.ax.plot(x, y, z, color='black', linewidth=0.3)
        elif argmax_index == 1:
            # 円盤生成
            geneNum = 200
            x = np.linspace(X_center-R, X_center+R, geneNum)
            z = np.linspace(Z_center-R, Z_center+R, geneNum)
            X, Z = np.meshgrid(x, z)
            if normalV[1] == 0:
                Y = Y_center - (normalV[0]*(X-X_center) +
                                normalV[2]*(Z-Z_center)) / 0.01
            else:
                Y = Y_center - (normalV[0]*(X-X_center) +
                                normalV[2]*(Z-Z_center)) / normalV[1]
            for i in range(geneNum):
                for j in range(geneNum):
                    if (X[i][j]-X_center)**2 + (Y[i][j]-Y_center)**2 + (Z[i][j]-Z_center)**2 > R**2:
                        Y[i][j] = np.nan
            #self.ax.quiver(X_center, Y_center, Z_center, normalV[0], normalV[1], normalV[2], color='black', length=50)
            self.ax.plot_wireframe(X, Y, Z, color='lightcyan', linewidth=0.3)
            theta = np.linspace(0, 2*np.pi, 100)
            x = R*np.cos(theta) + X_center
            y = np.zeros_like(theta) + Y_center
            z = R*np.sin(theta) + Z_center
            self.ax.plot(x, y, z, color='k', linewidth=0.3)
        elif argmax_index == 2:
            # 円盤生成
            geneNum = 200
            x = np.linspace(X_center-R, X_center+R, geneNum)
            y = np.linspace(Y_center-R, Y_center+R, geneNum)
            X, Y = np.meshgrid(x, y)
            if normalV[2] == 0:
                Z = Z_center - (normalV[0]*(X-X_center) +
                                normalV[1]*(Y-Y_center)) / 0.01
            else:
                Z = Z_center - (normalV[0]*(X-X_center) +
                                normalV[1]*(Y-Y_center)) / normalV[2]
            for i in range(geneNum):
                for j in range(geneNum):
                    if (X[i][j]-X_center)**2 + (Y[i][j]-Y_center)**2 + (Z[i][j]-Z_center)**2 > R**2:
                        Z[i][j] = np.nan
            #self.ax.quiver(X_center, Y_center, Z_center, normalV[0], normalV[1], normalV[2], color='black', length=50)
            self.ax.plot_wireframe(X, Y, Z, color='lightcyan', linewidth=0.3)
            theta = np.linspace(0, 2*np.pi, 100)
            x = R*np.cos(theta) + X_center
            y = R*np.sin(theta) + Y_center
            z = np.zeros_like(theta) + Z_center
            self.ax.plot(x, y, z, color='k', linewidth=0.3)

    # レンズ描画
    def plot_lens(self, params):
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
                self.ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)
            elif params[3] > 0:
                Xs1 = (params[3]**2-Ys1**2-Zs1**2)**0.5 - params[3]
                self.ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)
        elif argmax_index == 1:
            Xs = np.outer(np.sin(theta), np.sin(phi))
            Zs = np.outer(np.ones(np.size(theta)), np.cos(phi))
            Xs1 = params[2] * Xs
            Zs1 = params[2] * Zs
            if params[3] < 0:
                Ys1 = -(params[3]**2-Xs1**2-Zs1**2)**0.5 - params[3]
                self.ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)
            elif params[3] > 0:
                Ys1 = (params[3]**2-Xs1**2-Zs1**2)**0.5 - params[3]
                self.ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)
        elif argmax_index == 2:
            Xs = np.outer(np.sin(theta), np.sin(phi))
            Ys = np.outer(np.ones(np.size(theta)), np.cos(phi))
            Xs1 = params[2] * Xs
            Ys1 = params[2] * Ys
            if params[3] < 0:
                Zs1 = -(params[3]**2-Xs1**2-Ys1**2)**0.5 - params[3]
                self.ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)
            elif params[3] > 0:
                Zs1 = (params[3]**2-Xs1**2-Ys1**2)**0.5 - params[3]
                self.ax.plot_wireframe(Xs1+params[0][0], Ys1+params[0][1], Zs1+params[0][2], linewidth=0.1)

    # 放物線描画
    def plot_parabola(self, params):
        theta = np.linspace(0, 2*np.pi, 100)
        R = params[2]
        a = abs(1/(2*params[3]))
        x1 = R*np.cos(theta)
        z1 = R*np.sin(theta)
        X1, Z1 = np.meshgrid(x1, z1)
        Y1 = a*X1**2 + a*Z1**2
        for i in range(100):
            for j in range(100):
                if (X1[i][j])**2 + Z1[i][j]**2 > R**2:
                    Y1[i][j] = np.nan
                else:
                    Y1[i][j] = a*X1[i][j]**2 + a*Z1[i][j]**2
        self.ax.plot_wireframe(X1+params[0][0], Y1+params[0]
                        [1], Z1+params[0][2], color='b', linewidth=0.1)


    ## 光線計算
    # 受け取ったx,y,z座標から(x,y,z)の組を作る関数
    def makePoints(self, point0, point1, point2, shape0, shape1):
        result = [None]*(len(point0)+len(point1)+len(point2))
        result[::3] = point0
        result[1::3] = point1
        result[2::3] = point2
        result = np.array(result)
        result = result.reshape(shape0, shape1)
        return result

    # ベクトルの最大成分のインデックスを返す関数
    def max_index(self, vector):
        abs_vector = np.abs(vector)
        argmax_index = np.argmax(abs_vector)
        return argmax_index
    
    # ベクトルの最小成分のインデックスを返す関数
    def min_index(self, vector):
        abs_vector = np.abs(vector)
        argmin_index = np.argmin(abs_vector)
        return argmin_index

    # 光学素子の面情報を登録する関数
    def set_surface(self, surface):
        self.surface_pos = np.array(surface[0])  # [x, y, z]
        self.normalV_optical_element = np.array(surface[1])/np.linalg.norm(surface[1])  # [nV_x, nV_y, nV_z]
        self.normalV_refract_or_reflect = np.array(surface[1])/np.linalg.norm(surface[1])  # [nV_x, nV_y, nV_z]
        self.lens_or_parabola_R = surface[3]  # 曲率半径または焦点距離

    # 平板のレイトレーシング
    def raytrace_plane(self):
        centerV = self.surface_pos
        normalV = self.normalV_optical_element
        length_ray_start_dir = len(self.ray_start_dir)
        #print("length_ray_start_dir", length_ray_start_dir)
        if length_ray_start_dir == 3:
            nV = np.array(normalV)/np.linalg.norm(normalV)
            T = (np.dot(centerV, nV)-np.dot(self.ray_start_pos, nV)) / (np.dot(self.ray_start_dir, nV))
            self.ray_end_pos = self.ray_start_pos + T*self.ray_start_dir
            #print("VF:shape(ray_end_pos)!!!!!!!!!!!!!!!!!!!!!2", np.shape(self.ray_end_pos))
            self.optical_path += T*self.refractive_index_after
            self.normalV_refract_or_reflect = nV
        else:  # 光線群の場合
            nV = np.array(normalV)/np.linalg.norm(normalV)
            T = []
            for i in range(length_ray_start_dir):
                if np.dot(self.ray_start_dir[i], np.array([1,1,1])) == 0:
                    T_tmp = 0
                    print("VF, raytrace_plane, T_tmp", T_tmp)
                    T.append(T_tmp)
                else:
                    T_tmp = (np.dot(centerV, nV)-np.dot(np.array(self.ray_start_pos[i]), nV)) / (np.dot(np.array(self.ray_start_dir[i]), nV))
                    T.append(T_tmp)
            T = np.array(T)
            #T = [(np.dot(centerV, nV)-np.dot(self.ray_start_pos[i], nV)) / (np.dot(self.ray_start_dir[i], nV)) for i in range(length_ray_start_dir)]
            #print("VF:shape(ray_end_pos)!!!!!!!!!!!!!!!!!!!!!1", np.shape(self.ray_end_pos), np.shape(T), np.shape(self.ray_end_dir))
            self.ray_end_pos = self.ray_start_pos + np.array([V*T for V, T in zip(self.ray_start_dir, T)])
            #print("VF:shape(ray_end_pos)!!!!!!!!!!!!!!!!!!!!!2", np.shape(self.ray_end_pos), np.shape(T), np.shape(self.ray_end_dir), np.shape([V*T for V, T in zip(self.ray_start_dir, T)]))
            self.optical_path += np.array(T)*self.refractive_index_after
            self.normalV_refract_or_reflect = [nV]*length_ray_start_dir

    # 球面のレイトレーシング
    def raytrace_sphere(self):
        lens_pos = self.surface_pos
        lens_R = self.lens_or_parabola_R
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            tmp_V = np.zeros_like(self.ray_start_pos)
            tmp_index = self.max_index(self.ray_start_dir)
            tmp_V[tmp_index] = lens_R
            test_dot = np.dot(tmp_V, self.ray_start_dir)
            tmp_V = np.zeros_like(self.ray_start_pos)
            tmp_index = self.max_index(self.ray_start_dir)
            tmp_V[tmp_index] = lens_R
            shiftV = lens_pos - tmp_V
            #print("VFtest!!!!, shiftV =", shiftV)
            if test_dot > 0:  # 凹レンズ
                ray_pos = self.ray_start_pos - shiftV
                A = np.dot(self.ray_start_dir, self.ray_start_dir)
                B = np.dot(self.ray_start_dir, ray_pos)
                C = np.dot(ray_pos, ray_pos) - abs(lens_R)**2
                T = (-B + np.sqrt(B**2 - A*C)) / A
            elif test_dot < 0:  # 凸レンズ
                ray_pos = self.ray_start_pos - shiftV
                A = np.dot(self.ray_start_dir, self.ray_start_dir)
                B = np.dot(self.ray_start_dir, ray_pos)
                C = np.dot(ray_pos, ray_pos) - abs(lens_R)**2
                T = (-B - np.sqrt(B**2 - A*C)) / A
            else:
                T = 0
            self.ray_end_pos = self.ray_start_pos + T*self.ray_start_dir
            self.optical_path += np.array(T)*self.refractive_index_after
            self.normalV_refract_or_reflect = self.calcNormalV_sphere()
        else:  # 光線群の場合
            tmp_V = np.zeros_like(self.ray_start_pos)
            tmp_index = self.max_index(self.normalV_optical_element)
            tmp_V[:,tmp_index] = lens_R
            test_dot = np.dot(tmp_V[0], self.ray_start_dir[0])
            tmp_V = np.zeros(3)
            tmp_index = self.max_index(self.ray_start_dir[0])
            tmp_V[tmp_index] = lens_R
            shiftV = lens_pos - tmp_V
            #print("VFtest!!!!, shiftV =", shiftV)
            if test_dot > 0:  # 凹レンズ
                ray_pos = self.ray_start_pos - np.array([shiftV]*length_ray_start_dir)
                A = np.diag(np.dot(self.ray_start_dir, np.array(self.ray_start_dir).T))
                B = np.diag(np.dot(self.ray_start_dir, ray_pos.T))
                C = np.diag(np.dot(ray_pos, ray_pos.T)) - abs(lens_R)**2
                T = []
                for i in range(length_ray_start_dir):
                    if np.dot(self.ray_start_dir[i], np.array([1,1,1]))==0:
                        T_tmp = 0
                        T.append(T_tmp)
                    else:
                        T_tmp = (-B[i] + np.sqrt(B[i]**2 - A[i]*C[i])) / A[i]
                        T.append(T_tmp)
                T = np.array(T)
            elif test_dot < 0:  # 凸レンズ
                ray_pos = self.ray_start_pos - np.array([shiftV]*length_ray_start_dir)
                A = np.diag(np.dot(self.ray_start_dir, np.array(self.ray_start_dir).T))
                B = np.diag(np.dot(self.ray_start_dir, ray_pos.T))
                C = np.diag(np.dot(ray_pos, ray_pos.T)) - abs(lens_R)**2
                T = (-B - np.sqrt(B**2 - A*C)) / A
            else:
                T = np.zeros(length_ray_start_dir)
            self.ray_end_pos = self.ray_start_pos + [V*T for V, T in zip(self.ray_start_dir, T)]
            self.optical_path += np.array(T)*self.refractive_index_after
            self.normalV_refract_or_reflect = self.calcNormalV_sphere()

    # 球面の法線ベクトルを計算する関数
    def calcNormalV_sphere(self):
        surface_pos = self.surface_pos
        lens_R = self.lens_or_parabola_R
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            tmp_V = np.zeros(3)
            tmp_index = self.max_index(self.ray_start_dir)
            tmp_V[tmp_index] = lens_R
            normalV = self.ray_end_pos - surface_pos + tmp_V
            normalV = normalV / np.linalg.norm(normalV)
            return np.array(normalV)
        else:  # 光線群の場合
            tmp_V = np.zeros(3)
            tmp_index = self.max_index(self.ray_start_dir[0])
            tmp_V[tmp_index] = lens_R
            normalV = []
            for i in range(length_ray_start_dir):
                tmp_normalV = self.ray_end_pos[i] - surface_pos + tmp_V
                normalV.append(tmp_normalV/np.linalg.norm(tmp_normalV))
            return np.array(normalV)

    # 放物線のレイトレーシング
    def raytrace_parabola(self):
        parabola_pos = self.surface_pos
        parabola_R = self.lens_or_parabola_R
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            print("raytrace_parabola: インデックス自動化、未実装に注意")
            ray_pos = self.ray_start_pos - parabola_pos
            ray_dir = self.ray_start_dir
            if ray_dir[0]==0 and ray_dir[2]==0:
                a = 1/(2*parabola_R)
                T = a*(ray_pos[0]**2 - ray_pos[1]/a + ray_pos[0]**2) / ray_dir[1]
                self.ray_end_pos = self.ray_start_pos + T*self.ray_start_dir
                self.optical_path += np.array(T)*self.refractive_index_after
                self.normalV_refract_or_reflect = self.calcNormalV_parabola()
            else:
                if parabola_R<0:
                    a = -1/(2*parabola_R)
                    A = ray_dir[0]**2 + ray_dir[2]**2
                    B = ray_pos[0]*ray_dir[0] + ray_pos[2]*ray_dir[2] - ray_dir[1]/(2*a)
                    C = ray_pos[0]**2 + ray_pos[2]**2 - ray_pos[1]/a
                    T = (-B - np.sqrt(B**2 - A*C)) / A
                    self.ray_end_pos = self.ray_start_pos + T*self.ray_start_dir
                    self.optical_path += np.array(T)*self.refractive_index_after
                    self.normalV_refract_or_reflect = self.calcNormalV_parabola()
                else:
                    a = 1/(2*parabola_R)
                    A = ray_dir[0]**2 + ray_dir[2]**2
                    B = ray_pos[0]*ray_dir[0] + ray_pos[2]*ray_dir[2] - ray_dir[1]/(2*a)
                    C = ray_pos[0]**2 + ray_pos[2]**2 - ray_pos[1]/a
                    T = (-B + np.sqrt(B**2 - A*C)) / A
                    self.ray_end_pos = self.ray_start_pos + T*self.ray_start_dir
                    self.optical_path += np.array(T)*self.refractive_index_after
                    self.normalV_refract_or_reflect = self.calcNormalV_parabola()
        else:  # 光線群の場合
            print("raytrace_parabola: インデックス自動化、未実装に注意")
            ray_pos = self.ray_start_pos - parabola_pos
            ray_dir = self.ray_start_dir
            if ray_dir[0][0]==0 and ray_dir[0][2]==0:
                a = 1/(2*parabola_R)
                T = [a*(ray_pos[i][0]**2 - ray_pos[i][1]/a + ray_pos[i][0]**2) / ray_dir[i][1] for i in range(length_ray_start_dir)]
                self.ray_end_pos = self.ray_start_pos + T*self.ray_start_dir
                self.optical_path += np.array(T)*self.refractive_index_after
                self.normalV_refract_or_reflect = self.calcNormalV_parabola()
            else:
                if parabola_R<0:
                    a = -1/(2*parabola_R)
                    A = [ray_dir[i][0]**2 + ray_dir[i][2]**2 for i in range(length_ray_start_dir)]
                    B = [ray_pos[i][0]*ray_dir[i][0] + ray_pos[i][2]*ray_dir[i][2] - ray_dir[i][1]/(2*a) for i in range(length_ray_start_dir)]
                    C = [ray_pos[i][0]**2 + ray_pos[i][2]**2 - ray_pos[i][1]/a for i in range(length_ray_start_dir)]
                    T = [(-B[i] - np.sqrt(B[i]**2 - A[i]*C[i])) / A[i] for i in range(length_ray_start_dir)]
                    self.ray_end_pos = [self.ray_start_pos[i] + T[i]*self.ray_start_dir[i] for i in range(length_ray_start_dir)]
                    self.optical_path += [np.array(T[i])*self.refractive_index_after for i in range(length_ray_start_dir)]
                    self.normalV_refract_or_reflect = self.calcNormalV_parabola()
                else:
                    a = 1/(2*parabola_R)
                    A = [ray_dir[i][0]**2 + ray_dir[i][2]**2 for i in range(length_ray_start_dir)]
                    B = [ray_pos[i][0]*ray_dir[i][0] + ray_pos[i][2]*ray_dir[i][2] - ray_dir[i][1]/(2*a) for i in range(length_ray_start_dir)]
                    C = [ray_pos[i][0]**2 + ray_pos[i][2]**2 - ray_pos[i][1]/a for i in range(length_ray_start_dir)]
                    T = [(-B[i] + np.sqrt(B[i]**2 - A[i]*C[i])) / A[i] for i in range(length_ray_start_dir)]
                    self.ray_end_pos = [self.ray_start_pos[i] + T[i]*self.ray_start_dir[i] for i in range(length_ray_start_dir)]
                    self.optical_path += [np.array(T[i])*self.refractive_index_after for i in range(length_ray_start_dir)]
                    self.normalV_refract_or_reflect = self.calcNormalV_parabola()

    # 放物線の法線ベクトルを計算する関数
    def calcNormalV_parabola(self):
        parabola_pos = self.surface_pos
        parabola_R = self.lens_or_parabola_R
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            tmp_index = self.max_index(self.normalV_optical_element)  # 方向の計算に使う
            a = abs(1/(2*parabola_R))
            ray_pos = self.ray_end_pos
            normalV = 2*a*(ray_pos - parabola_pos)
            normalV[tmp_index] = -1
            normalV = normalV/np.linalg.norm(normalV)
            return normalV
        else:  # 光線群の場合
            tmp_index = self.max_index(self.normalV_optical_element)  # 方向の計算に使う
            a = abs(1/(2*parabola_R))
            ray_pos = self.ray_end_pos
            normalV = 2*a*(ray_pos - parabola_pos)
            normalV[:,tmp_index] = -1
            normalV = np.array([V/np.linalg.norm(V) for V in normalV])
            return normalV

    # 反射計算
    def reflect(self):
        normalV = self.normalV_refract_or_reflect
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            ray_dir = np.array(self.ray_start_dir)/np.linalg.norm(self.ray_start_dir)
            normalV = np.array(normalV)/np.linalg.norm(normalV)
            outRayV = ray_dir - 2*(np.dot(ray_dir, normalV))*normalV
            # 正規化
            outRayV = outRayV/np.linalg.norm(outRayV)
            self.ray_end_dir = outRayV
        else:  # 光線群の場合
            ray_dir = [V/np.linalg.norm(V) for V in self.ray_start_dir]
            normalV = [V/np.linalg.norm(V) for V in normalV]
            #dot_tmp = np.dot(ray_dir, normalV)
            dot_tmp = [np.dot(ray_dir[i], normalV[i]) for i in range(length_ray_start_dir)]
            dot_tmp2 = [np.dot(dot_tmp[i], normalV[i]) for i in range(length_ray_start_dir)]
            outRayV = np.array(ray_dir) - 2*np.array(dot_tmp2)
            # 正規化
            outRayV = [V/np.linalg.norm(V) for V in outRayV]
            self.ray_end_dir = outRayV

    # スネルの法則
    def refract(self):
        normalV = self.normalV_refract_or_reflect
        Nin = self.refractive_index_before
        Nout = self.refractive_index_after
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            # 正規化
            ray_dir = self.ray_start_dir/np.linalg.norm(self.ray_start_dir)
            normalV = normalV/np.linalg.norm(normalV)
            if np.dot(ray_dir, normalV) <= 0:
                #print("内積が負です")
                # 係数A
                A = Nin/Nout
                # 入射角
                cos_t_in = abs(np.dot(ray_dir, normalV))
                # 量子化誤差対策
                if cos_t_in < -1.:
                    cos_t_in = -1.
                elif cos_t_in > 1.:
                    cos_t_in = 1.
                # スネルの法則
                sin_t_in = np.sqrt(1.0 - cos_t_in**2)
                sin_t_out = sin_t_in*A
                if sin_t_out > 1.0:
                    # 全反射する場合
                    return np.zeros(3)
                cos_t_out = np.sqrt(1 - sin_t_out**2)
                # 係数B
                B = A*cos_t_in - cos_t_out
                # 出射光線の方向ベクトル
                outRayV = A*ray_dir + B*normalV
                # 正規化
                outRayV = outRayV/np.linalg.norm(outRayV)
            else:
                #print("内積が正です")
                # 係数A
                A = Nin/Nout
                # 入射角
                cos_t_in = abs(np.dot(ray_dir,normalV))
                #量子化誤差対策
                if cos_t_in<-1.:
                    cos_t_in = -1.
                elif cos_t_in>1.:
                    cos_t_in = 1.
                # スネルの法則
                sin_t_in = np.sqrt(1.0 - cos_t_in**2)
                sin_t_out = sin_t_in*A
                if sin_t_out>1.0:
                    #全反射する場合
                    return np.zeros(3)
                cos_t_out = np.sqrt(1 - sin_t_out**2)
                # 係数B
                B = -A*cos_t_in + cos_t_out
                # 出射光線の方向ベクトル
                outRayV = A*ray_dir + B*normalV
                # 正規化
                outRayV = outRayV/np.linalg.norm(outRayV)
            self.ray_end_dir = outRayV
        else:
            # 正規化
            ray_dir = [V/np.linalg.norm(V) for V in self.ray_start_dir]
            normalV = [V/np.linalg.norm(V) for V in normalV]
            tmp_V = np.zeros(3)
            if len(normalV) == 3:
                tmp_index = self.max_index(normalV)
            else:
                tmp_index = self.max_index(normalV[0])
            tmp_V[tmp_index] = 1.
            test_dot = np.dot(ray_dir[0], tmp_V)
            #if test_dot <= 0:
            if np.dot(ray_dir[0], normalV[0]) <= 0:
                #print("内積が負です")
                # 係数A
                A = Nin/Nout
                # 入射角
                cos_t_in = np.abs(np.diag(np.dot(ray_dir,np.array(normalV).T)))
                outRayV = []
                for i in range(length_ray_start_dir):
                    i_cos_t_in = cos_t_in[i]
                    i_ray_dir = ray_dir[i]
                    i_normalV = normalV[i]
                    # 量子化誤差対策
                    if i_cos_t_in < -1.:
                        i_cos_t_in = -1.
                    elif i_cos_t_in > 1.:
                        i_cos_t_in = 1.
                    # スネルの法則
                    sin_t_in = np.sqrt(1.0 - i_cos_t_in**2)
                    sin_t_out = sin_t_in*A
                    if sin_t_out > 1.0:
                        # 全反射する場合
                        return np.zeros(3)
                    cos_t_out = np.sqrt(1 - sin_t_out**2)
                    # 係数B
                    B = A*i_cos_t_in - cos_t_out
                    # 出射光線の方向ベクトル
                    tmp_outRayV = A*i_ray_dir + B*i_normalV
                    # 正規化
                    tmp_outRayV = tmp_outRayV/np.linalg.norm(tmp_outRayV)
                    outRayV.append(tmp_outRayV)
            else:
                #print("内積が正です")
                # 係数A
                A = Nin/Nout
                # 入射角
                cos_t_in = np.abs(np.diag(np.dot(ray_dir,np.array(normalV).T)))
                outRayV = []
                for i in range(length_ray_start_dir):
                    i_cos_t_in = cos_t_in[i]
                    i_ray_dir = ray_dir[i]
                    i_normalV = normalV[i]
                    #量子化誤差対策
                    if i_cos_t_in<-1.:
                        i_cos_t_in = -1.
                    elif i_cos_t_in>1.:
                        i_cos_t_in = 1.
                    # スネルの法則
                    sin_t_in = np.sqrt(1.0 - i_cos_t_in**2)
                    sin_t_out = sin_t_in*A
                    if sin_t_out>1.0:
                        #全反射する場合
                        return np.zeros(3)
                    cos_t_out = np.sqrt(1 - sin_t_out**2)
                    # 係数B
                    B = -A*i_cos_t_in + cos_t_out
                    # 出射光線の方向ベクトル
                    tmp_outRayV = A*i_ray_dir + B*i_normalV
                    # 正規化
                    tmp_outRayV = tmp_outRayV/np.linalg.norm(tmp_outRayV)
                    outRayV.append(tmp_outRayV)
            self.ray_end_dir = outRayV


    # 焦点距離を計算する関数
    def calcFocalLength(self, ray_start_pos_init):
        length_ray_dir = len(self.ray_end_pos)
        if length_ray_dir == 3:
            argmin_index = self.min_index(self.ray_end_dir)
            tmp_V = -1.*self.ray_end_dir*ray_start_pos_init[argmin_index]/self.ray_end_dir[argmin_index]
            argmax_index = self.max_index(tmp_V)
            focal_length = tmp_V[argmax_index]
            return focal_length
        else:  # 光線が複数の場合
            argmin_index = self.min_index(self.ray_end_dir[0])
            focal_length = []
            for i in range(length_ray_dir):
                tmp_V = -1.*self.ray_end_dir[i]*ray_start_pos_init[i][argmin_index]/self.ray_end_dir[i][argmin_index]
                argmax_index = self.max_index(tmp_V)
                focal_length.append(tmp_V[argmax_index])
            # 並び替え
            focal_length.sort()
            return focal_length[0]

    # 焦点位置を計算する関数
    def calcFocalPos(self, ray_start_pos_init):
        length_ray_dir = len(self.ray_end_pos)
        if length_ray_dir == 3:
            argmin_index = self.min_index(self.ray_end_dir)
            tmp_V = -1.*self.ray_end_dir*self.ray_end_pos[argmin_index]/self.ray_end_dir[argmin_index]
            focal_point = tmp_V + self.ray_end_pos
            argmax_index = self.max_index(tmp_V)
            print("!!!!正確な焦点位置を得るには近軸光線を計算する必要があります!!!!")
            return focal_point[argmax_index]
        else:
            argmin_index = self.min_index(self.ray_end_dir[0])
            focal_point_list = []
            for i in range(length_ray_dir):
                tmp_V = -1.*self.ray_end_dir[i]*self.ray_end_pos[i][argmin_index]/self.ray_end_dir[i][argmin_index]
                focal_point_list.append(tmp_V + self.ray_end_pos[i])
            r_list = []
            focus_pos = []
            for i in range(length_ray_dir):
                r_list.append(np.sqrt(ray_start_pos_init[i][1]**2 + ray_start_pos_init[i][2]**2))
                focus_pos.append(focal_point_list[i][0])
            # 並び替え
            r_list, focus_pos = zip(*sorted(zip(r_list, focus_pos)))
            #print("focus_pos = ", focus_pos)
            focal_point = focus_pos[0]
            print("!!!!正確な焦点位置を得るには近軸光線を計算する必要があります!!!!")
            return focal_point


    # ２点の位置ベクトルから直線を引く関数
    def plotLineBlue(self):
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            startPointV = self.ray_start_pos
            endPointV = self.ray_end_pos
            startX = startPointV[0]
            startY = startPointV[1]
            startZ = startPointV[2]
            endX = endPointV[0]
            endY = endPointV[1]
            endZ = endPointV[2]
            self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                    'o-', ms='2', linewidth=0.5, color='blue')
        else:
            for i in range(length_ray_start_dir):
                startPointV = self.ray_start_pos[i]
                endPointV = self.ray_end_pos[i]
                startX = startPointV[0]
                startY = startPointV[1]
                startZ = startPointV[2]
                endX = endPointV[0]
                endY = endPointV[1]
                endZ = endPointV[2]
                self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                        'o-', ms='2', linewidth=0.5, color='blue')

    def plotLineGreen(self):
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            startPointV = self.ray_start_pos
            endPointV = self.ray_end_pos
            startX = startPointV[0]
            startY = startPointV[1]
            startZ = startPointV[2]
            endX = endPointV[0]
            endY = endPointV[1]
            endZ = endPointV[2]
            self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                    'o-', ms='2', linewidth=0.5, color='green')
        else:
            for i in range(length_ray_start_dir):
                startPointV = self.ray_start_pos[i]
                endPointV = self.ray_end_pos[i]
                startX = startPointV[0]
                startY = startPointV[1]
                startZ = startPointV[2]
                endX = endPointV[0]
                endY = endPointV[1]
                endZ = endPointV[2]
                self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                        'o-', ms='2', linewidth=0.5, color='green')

    def plotLineRed(self):
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            startPointV = self.ray_start_pos
            endPointV = self.ray_end_pos
            startX = startPointV[0]
            startY = startPointV[1]
            startZ = startPointV[2]
            endX = endPointV[0]
            endY = endPointV[1]
            endZ = endPointV[2]
            self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                    'o-', ms='2', linewidth=0.5, color='r')
        else:
            for i in range(length_ray_start_dir):
                startPointV = self.ray_start_pos[i]
                endPointV = self.ray_end_pos[i]
                startX = startPointV[0]
                startY = startPointV[1]
                startZ = startPointV[2]
                endX = endPointV[0]
                endY = endPointV[1]
                endZ = endPointV[2]
                self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                        'o-', ms='2', linewidth=0.5, color='r')

    def plotLineOrange(self):
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            startPointV = self.ray_start_pos
            endPointV = self.ray_end_pos
            startX = startPointV[0]
            startY = startPointV[1]
            startZ = startPointV[2]
            endX = endPointV[0]
            endY = endPointV[1]
            endZ = endPointV[2]
            self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                    'o-', ms='2', linewidth=0.5, color='orange')
        else:
            for i in range(length_ray_start_dir):
                startPointV = self.ray_start_pos[i]
                endPointV = self.ray_end_pos[i]
                startX = startPointV[0]
                startY = startPointV[1]
                startZ = startPointV[2]
                endX = endPointV[0]
                endY = endPointV[1]
                endZ = endPointV[2]
                self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                        'o-', ms='2', linewidth=0.5, color='orange')
 
    def plotFourBeamLine(self, i):
        if i == 0:
            self.plotLineBlue()
        elif i == 1:
            self.plotLineGreen()
        elif i == 2:
            self.plotLineRed()
        elif i == 3:
            self.plotLineOrange()
        else:
            print("fourBeamPlotLine, iの値が不正です")

    def plotLinePurple(self, startPointV, endPointV):
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            startX = startPointV[0]
            startY = startPointV[1]
            startZ = startPointV[2]
            endX = endPointV[0]
            endY = endPointV[1]
            endZ = endPointV[2]
            self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                    'o-', ms='2', linewidth=0.5, color='purple')
        else:
            for i in range(length_ray_start_dir):
                startPointV = self.ray_start_pos[i]
                endPointV = self.ray_end_pos[i]
                startX = startPointV[0]
                startY = startPointV[1]
                startZ = startPointV[2]
                endX = endPointV[0]
                endY = endPointV[1]
                endZ = endPointV[2]
                self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                        'o-', ms='2', linewidth=0.5, color='purple')

    def plotLineBlack(self, startPointV, endPointV):
        length_ray_start_dir = len(self.ray_start_pos)
        if length_ray_start_dir == 3:
            startX = startPointV[0]
            startY = startPointV[1]
            startZ = startPointV[2]
            endX = endPointV[0]
            endY = endPointV[1]
            endZ = endPointV[2]
            self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                    'o-', ms='2', linewidth=0.5, color='black')
        else:
            for i in range(length_ray_start_dir):
                startPointV = self.ray_start_pos[i]
                endPointV = self.ray_end_pos[i]
                startX = startPointV[0]
                startY = startPointV[1]
                startZ = startPointV[2]
                endX = endPointV[0]
                endY = endPointV[1]
                endZ = endPointV[2]
                self.ax.plot([startX, endX], [startY, endY], [startZ, endZ],
                        'o-', ms='2', linewidth=0.5, color='black')