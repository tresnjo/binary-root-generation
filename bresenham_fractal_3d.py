import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def bresenham3D(x1, x2, y1, y2, z1, z2, th, matrix):

    ''' ================= VARIABLE DEFINITIONS =================
                x1, y1 - the initial coordinates
                x2, y2 - the final coordinates
                th - the thickness of the line
                matrix - matrix for coordinates                '''

    # Rounding all coordinates to nearest integer 
    
    x2 = int(round(x2))
    x1 = int(round(x1))
    y2 = int(round(y2))
    y1 = int(round(y1))
    z1 = int(round(z1))
    z2 = int(round(z2))

    # Calculating the differences between the two coordinates
    xdif = x2-x1
    ydif = y2-y1
    zdif = z2-z1
    
    # Absolute values of differences
    dx = abs(x2-x1)
    dy = abs(y2-y1)
    dz = abs(z2-z1)

    # Setting signs of increments in x, y and z - direction
    xs = 1 if xdif > 0 else -1
    ys = 1 if ydif > 0 else -1
    zs = 1 if zdif > 0 else -1

    dx_2 = 2*dx
    dy_2 = 2*dy
    dz_2 = 2*dz

    min_thick = 1

    # lägg till case dx = 0 dy = 0 men dz != 0 etc
    # Case I: dx > dy and dx > dz
    if (dx >= dy and dx >= dz): 
            p1 = dy_2 - dx
            p2 = dz_2- dx
            while (x1 != x2):

                if th == 1:
                    matrix[y1][x1][z1] = 1
                else:
                    for y_th in range(0, 2*th+1):
                        max_t_z = max(int(np.sqrt(th**2-(y_th-th)**2)), min_thick)
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
                if th == 1:
                    matrix[y1][x1][z1] = 1
                else:
                    for x_th in range(0, 2*th+1):
                        max_t_z = max(int(np.sqrt(th**2-(x_th-th)**2)), min_thick)
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
            if th == 1:
                 matrix[y1][x1][z1] = 1
            else:
                for x_th in range(0, 2*th+1):
                    max_t_y = max(int(np.sqrt(th**2-(x_th-th)**2)), min_thick)
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

    # non uniformity? removing branches after shuffling?
    # adding a minimum thickness?
    # FIX RETURN ERROR IF INDEX OUT OF BOUNDS

def FraktalT(n, r, phi, chi, xb, yb, zb, t, tr, dim, pixels):

    ''' ================= VARIABLE DEFINITIONS =================
                        n - number of iterations

             r - scale factor (different scale-factors take                    
            place, r may be vector in this case the lengths 
                    of r and phi must be equal)

        phi - vector of angles in fractal generator (calculated 
        relative to the horizontal axis connecting end points 
                            of generator)
        chi - vector of polar angles in fractal generator (calculated 
        relative to the horizontal axis connecting end points 
                            of generator)
    
                 xb,yb and zb - start and end coordinates of trunk 
                        t - radius of root in mm

        tr - scale factor for thickness (different scale factors 
             can take place BUT lengths must mach with 
                         lengths of r and phi)

            dim - dimensions of domain in mm 
            pixels - number of pixels along x and y axis        '''

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
    for i in range(1, n):
        z = 0
        for j in range(0, mN**(i - 1)):
            for k in range(0, mN):
                for m in range(0, mN**(n - i)):

                    a = np.sqrt((A[i - 1, z, 0] - A[i, z, 0])**2 + (A[i - 1, z, 1] - A[i, z, 1])**2 + (A[i - 1, z, 2] - A[i, z, 2])**2)
                    b = np.sqrt((A[i, z, 0] - A[i - 1, z, 0])**2 + (A[i, z, 1] - A[i - 1, z, 1])**2)
                    theta1 = np.arccos((A[i, z, 2] - A[i - 1, z, 2]) / a)

                    if A[i, z, 0] == A[i - 1, z, 0] and A[i, z, 1] == A[i - 1, z, 1]:
                        k2 = 0
                        k1 = 1
                    else:
                        k2 = (A[i, z, 1] - A[i - 1, z, 1]) / b
                        k1 = (A[i, z, 0] - A[i - 1, z, 0]) / b
                        
                    B = np.array([[k1 * np.cos(theta1), -k2, np.sin(theta1) * k1],
                                  [k2 * np.cos(theta1), k1, k2 * np.sin(theta1)],
                                  [-np.sin(theta1), 0, np.cos(theta1)]])

                    x2 = a * rM[k] * np.sin(chi[k]) * np.cos(phi[k])
                    y2 = a * rM[k] * np.sin(chi[k]) * np.sin(phi[k])
                    z2 = a * rM[k] * np.cos(chi[k])

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
                matrix = bresenham3D(A[j-1, i-1, 0], A[j, i-1, 0], A[j-1, i-1, 1], A[j, i-1, 1], A[j-1, i-1,2], A[j, i-1,2] ,int(t*tr**(i-1)), matrix)
            z += 1

    # Visualizing end branches
    for i in range(mN**n):
        matrix = bresenham3D(A[n-1, i-1, 0], A[n, i-1, 0], A[n-1, i-1, 1], A[n , i-1, 1],A[n-1, i-1,2], A[n, i-1,2], int(t*tr**(n)), matrix)

     # Visualizing the trunk
    matrix = bresenham3D(xb[0]/delta_x + pixels[1]//2, xb[1]/delta_x + pixels[1]//2, yb[0]/delta_y + pixels[0]//2, yb[1]/delta_y + pixels[0]//2, zb[0]/delta_z, zb[1]/delta_z,int(t), matrix)
    return matrix

matrix = FraktalT(3,[0.8, 0.6, 0.8, 0.6],[np.pi/4, -np.pi/4, 3*np.pi/4, 5*np.pi/4, ],[np.pi/3, np.pi/6,np.pi/6, np.pi/3], [0,0], [0,0], [0,100], 5, 0.8,[500,500,500], [120,120,120])

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection="3d")
ax.invert_xaxis()
ax.voxels(matrix, facecolor = 'black', edgecolor = 'white',shade=None)
plt.show()

