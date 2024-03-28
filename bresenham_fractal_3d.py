import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random as random
import matplotlib.colors


### Configuration Options ##
DEPTH = 1                                                         # Fractal depth
BRANCHING_RATIO = [0.8, 0.6, 0.8, 0.6]                              # Length ratio bw parent and daugther branch
BRANCHING_POLAR = [np.pi/4, -np.pi/4, 3*np.pi/4, 5*np.pi/4, ]       # Polar angles
BRANCHING_AZIMUTHAL = [np.pi/3, np.pi/6, np.pi/6, np.pi/3]          # Azimuthal angles

ROOT_START = [0,0,0]                                                # Start of root
ROOT_END = [0,0,100]                                                # End of root
THICKNESS_ROOT = 4                                                  # Thickness of the root
MIN_THICKNESS_BRANCH = 2                                            # Minimum thickness of the branch

THICKNESS_RATIO = 0.8                                               # Thickness ratio bw parent and daugther branch     
DOMAIN_DIMENSIONS = [500, 500, 500]                                 # Dimensions of the domain 
DOMAIN_PIXELS = [120, 120, 120]                                     # Number of pixels in binary domain            

RANDOM_MODE = False                                                 # Boolean for random fractal
SEED = 412                                                          # Seed for random fractal generation
ANGLE_VARIATION_POLAR = np.pi/8                                     # Maximum deviation in polar angles
ANGLE_VARIATION_AZIMUTHAL = np.pi/8                                 # Maximum deviation in azimuthal angles
BRANCHING_VARIATION = 0.1                                           # Maximum deviation in branching ratio
THICKNESS_VARIATION = 2                                            # Maximum deviation in thickness

random.seed(SEED)

def bresenham3D(x1, x2, y1, y2, z1, z2, th, matrix):

    # Rounding all coordinates to nearest integer 
    x1, x2 =int(round(x1)), int(round(x2))
    y1, y2 =int(round(y1)), int(round(y2))
    z1, z2 =int(round(z1)), int(round(z2))

    # Calculating the differences between the two coordinates
    xdif, ydif, zdif = x2-x1, y2-y1, z2-z1
    
    # Absolute values of differences
    dx,dy, dz = abs(x2-x1), abs(y2-y1), abs(z2-z1)

    # Setting signs of increments in x, y and z - direction
    xs = 1 if xdif > 0 else -1
    ys = 1 if ydif > 0 else -1
    zs = 1 if zdif > 0 else -1

    dx_2,dy_2,dz_2 = 2*dx,2*dy,2*dz

    if RANDOM_MODE == True:
        th += random.uniform(-1,1)*THICKNESS_VARIATION
        th = int(th)

    # Case I: dx > dy and dx > dz
    if (dx >= dy and dx >= dz): 
            p1 = dy_2 - dx
            p2 = dz_2- dx
            while (x1 != x2):
                
                if y1 >= DOMAIN_PIXELS[1] or z1 >= DOMAIN_PIXELS[2] or x1 >= DOMAIN_PIXELS[0] or x1 < 0 or y1 < 0 or z1 < 0:
                        print("Error. Coordinate out of bounds.")
                        continue
                else:
                    if th == 1:
                        matrix[y1][x1][z1] = 1
                    else:
                        for y_th in range(0, 2*th+1):
                            max_t_z = max(int(np.sqrt(th**2-(y_th-th)**2)), MIN_THICKNESS_BRANCH)
                            for z_th in range(0,2*max_t_z+1):
                                matrix[y1-th+y_th][x1][z1-max_t_z+z_th] = 1

                
                if (p1 >= 0):
                    y1 += ys
                    p1 -= dx_2
                if (p2 >= 0):
                    z1 += zs
                    p2 -= dx_2
                    
                p1 += dy_2
                p2 += dz_2
                x1 += xs

    # Case II: dy > dx and dy > dz
    elif (dy >= dx and dy >= dz): 
            p1 = dx_2 - dy
            p2 = dz_2- dy
            while (y1 != y2):
                 
                if y1 >= DOMAIN_PIXELS[1] or z1 >= DOMAIN_PIXELS[2] or x1 >= DOMAIN_PIXELS[0] or x1 < 0 or y1 < 0 or z1 < 0:
                        print("Error. Coordinate out of bounds.")
                        continue
                else:
                    if th == 1:
                        matrix[y1][x1][z1] = 1
                    else:
                        for x_th in range(0, 2*th+1):
                            max_t_z = max(int(np.sqrt(th**2-(x_th-th)**2)), MIN_THICKNESS_BRANCH)
                            for z_th in range(0,2*max_t_z+1):
                                matrix[y1][x1-th+x_th][z1-max_t_z+z_th] = 1
                if (p1 >= 0):
                    x1 += xs
                    p1 -= dy_2
                if (p2 >= 0):
                    z1 += zs
                    p2 -= dy_2
                p1 += dx_2
                p2 += dz_2
                y1 += ys

    # Case III: dz > dx and dz > dy
    else:
        p1 = dy_2 - dz
        p2 = dx_2- dz
        while (z1 != z2):
            if y1 >= DOMAIN_PIXELS[1] or z1 >= DOMAIN_PIXELS[2] or x1 >= DOMAIN_PIXELS[0] or x1 < 0 or y1 < 0 or z1 < 0:
                        print("Error. Coordinate out of bounds.")
                        continue
            else:
                if th == 1:
                    matrix[y1][x1][z1] = 1
                else:
                    for x_th in range(0, 2*th+1):
                        max_t_y = max(int(np.sqrt(th**2-(x_th-th)**2)), MIN_THICKNESS_BRANCH)
                        for y_th in range(0,2*max_t_y+1):
                            matrix[y1-max_t_y + y_th][x1-th+x_th][z1] = 1
            if (p1 >= 0):
                y1 += ys
                p1 -= dz_2
            if (p2 >= 0):
                x1 += xs
                p2 -= dz_2
            p1 += dy_2
            p2 += dx_2
            z1 += zs

    return matrix

def FraktalT(n, r, phi, chi, xb, yb, zb, t, tr, dim, pixels, rnd_phi, rnd_chi, rnd_branch):

    # Calculating voxel size given pixel and domain size 
    delta_x = dim[0]/pixels[1]
    delta_y = dim[1]/pixels[0]
    delta_z = dim[2]/pixels[2]

    # Initializing matrix for domain
    matrix = np.zeros((pixels[1], pixels[0], pixels[2]))

    # Initializnig matrix for scale factors
    mN = len(phi)
    if len(r) == 1:
        rM = np.ones(mN) * r
    elif len(r) == mN:
        rM = np.array(r)
    else:
        print('The sizes of scale vector of lengths and vector of angles in fractal`s generator don`t equal')
        return
    
    # Initialize matrix for coordinates
    A = np.ones((n + 1, mN**n, 3))

    # Trunk coordinates
    A[0,:,0] = np.ones((1,mN**n))*xb[0]/delta_x + pixels[1]//2
    A[0,:,1] = np.ones((1,mN**n))*yb[0]/delta_y+ pixels[0]//2
    A[0,:,2] = np.ones((1,mN**n))*zb[0]/delta_z
    A[1,:,0] =  np.ones((1,mN**n))*xb[1]/delta_x+ pixels[1]//2
    A[1,:,1]  =  np.ones((1,mN**n))*yb[1]/delta_y+ pixels[0]//2
    A[1,:,2]  =  np.ones((1,mN**n))*zb[1]/delta_z

    NC = np.zeros((3, mN))
    if RANDOM_MODE == False:
        for i in range(1, n):
            z = 0
            for j in range(0, mN**(i - 1)):
                for k in range(0, mN):
                    for m in range(0, mN**(n - i)):

                        a = np.sqrt((A[i - 1, z, 0] - A[i, z, 0])**2 
                                    + (A[i - 1, z, 1] - A[i, z, 1])**2 +
                                    (A[i - 1, z, 2] - A[i, z, 2])**2)
    
                        x2 = a * rM[k] * np.sin(chi[k]) * np.cos(phi[k])
                        y2 = a * rM[k] * np.sin(chi[k]) * np.sin(phi[k])
                        z2 = a * rM[k] * np.cos(chi[k])

                        NC[:, k] = np.array([x2, y2, z2]).T + A[i, z].T
            
                        # Define following coordinates
                        A[i + 1, z,0] = NC[0,k] 
                        A[i + 1, z,1] = NC[1,k]
                        A[i + 1, z,2] = NC[2,k]

                        z += 1
    else:
         for i in range(1, n):
            z = 0
            for j in range(0, mN**(i - 1)):
                for k in range(0, mN):
                    for m in range(0, mN**(n - i)):

                        a = np.sqrt((A[i - 1, z, 0] - A[i, z, 0])**2 
                                    + (A[i - 1, z, 1] - A[i, z, 1])**2 +
                                    (A[i - 1, z, 2] - A[i, z, 2])**2)
                            
                        x2 = a * (rM[k]+np.random.uniform(-1,1)*rnd_branch) * np.sin(chi[k]+np.random.uniform(-1,1)*rnd_chi) * np.cos(phi[k]+np.random.uniform(-1,1)*rnd_phi)
                        y2 = a * (rM[k]+np.random.uniform(-1,1)*rnd_branch) * np.sin(chi[k]+np.random.uniform(-1,1)*rnd_chi) * np.sin(phi[k]+np.random.uniform(-1,1)*rnd_phi)
                        z2 = a * (rM[k]+np.random.uniform(-1,1)*rnd_branch) * np.cos(chi[k]+np.random.uniform(-1,1)*rnd_chi)

                        NC[:, k] = np.array([x2, y2, z2]).T + A[i, z].T
            
                        # Define following coordinates
                        A[i + 1, z,0] = NC[0,k] 
                        A[i + 1, z,1] = NC[1,k]
                        A[i + 1, z,2] = NC[2,k]

                        z += 1
    # Visualizing lines using Bresenhams line algorithm
    for i in range(1, mN**n+1, mN):
        z = 1
        for k in range(1,mN):
            for j in range(z, n ):
                matrix = bresenham3D(A[j-1, i-1, 0], A[j, i-1, 0], A[j-1, i-1, 1], 
                                     A[j, i-1, 1], A[j-1, i-1,2], A[j, i-1,2] ,
                                     int(round(t*tr**(i-1))), matrix)
            z += 1
    # Visualizing end branches
    for i in range(mN**n):
        matrix = bresenham3D(A[n-1, i-1, 0], A[n, i-1, 0], A[n-1, i-1, 1], 
                             A[n , i-1, 1],A[n-1, i-1,2], A[n, i-1,2], 
                             int(round(t*tr**(n))), matrix)

     # Visualizing the trunk
    matrix = bresenham3D(xb[0]/delta_x + pixels[1]//2, xb[1]/delta_x + pixels[1]//2, 
                         yb[0]/delta_y + pixels[0]//2, yb[1]/delta_y + pixels[0]//2, 
                         zb[0]/delta_z, zb[1]/delta_z,
                         int(round(t)), matrix)
    return matrix

matrix = FraktalT(DEPTH,BRANCHING_RATIO,BRANCHING_POLAR,BRANCHING_AZIMUTHAL, 
                  [ROOT_START[0], ROOT_END[0]], [ROOT_START[1], ROOT_END[1]],[ROOT_START[2], ROOT_END[2]], 
                  THICKNESS_ROOT, THICKNESS_RATIO, DOMAIN_DIMENSIONS, DOMAIN_PIXELS,
                  ANGLE_VARIATION_POLAR, ANGLE_VARIATION_AZIMUTHAL, BRANCHING_VARIATION)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection="3d")
ax.invert_xaxis()
ax.set_axis_off()
ax.voxels(matrix, facecolor = 'brown', edgecolor = 
          'black', shade=None)
plt.show()

