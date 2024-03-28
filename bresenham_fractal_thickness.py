import numpy as np
import matplotlib.pyplot as plt

def bresenham2D(x1, x2, y1, y2, th, matrix):

    ''' ================= VARIABLE DEFINITIONS =================
                x1, y1 - the initial coordinates
                x2, y2 - the final coordinates
                th - the thickness of the line
                matrix - matrix for coordinates      
                                                                '''

    # Minimum thickness of any branch
    min_thickness = 2
    th = max(th, min_thickness)

    # Rounding all coordinates to nearest integer 
    x2 = int(round(x2))
    x1 = int(round(x1))
    y2 = int(round(y2))
    y1 = int(round(y1))

    # Calculating the differences between the two coordinates
    xdif = x2-x1
    ydif = y2-y1
    
    # Absolute values of differences
    dx = abs(x2-x1)
    dy = abs(y2-y1)

    # Setting initial point as 1 in matrix and looping over thickness
    for i in range(0,2*th+1):
        if dx > dy:
             matrix[y1+i-th][x1] = 1
        elif dx < dy:
            matrix[y1][x1+ i - th] = 1
        else:
            matrix[y1+i-th][x1] = 1
            
    # Setting signs of increments in x, y - direction
    xs = 1 if xdif > 0 else -1
    ys = 1 if ydif > 0 else -1
        
    # Case I: dx > dy
    if (dx > dy): 
            p = 2 * dy - dx
            while (x1 != x2):
                x1 += xs
                if (p >= 0):
                    y1 += ys
                    p -= 2 * dx
                p += 2 * dy
                # looping over thickness
                for i in range(0, 2*th+1):
                    matrix[y1+i-th][x1] = 1

    # Case II: dx < dy
    else:
            p = 2 * dx-dy
            while(y1 != y2):
                y1 += ys
                if (p >= 0):
                    x1 += xs
                    p -= 2 * dy
                p += 2 * dx
                # looping over thickness
                for i in range(0, 2*th+1):
                    matrix[y1][x1+ i - th] = 1

    return matrix

    # non uniformity? removing branches after shuffling?
    # adding a minimum thickness?
    # FIX RETURN ERROR IF INDEX OUT OF BOUNDS

def FraktalT(n, r, phi, xb, yb, t, tr, dim, pixels):

    ''' ================= VARIABLE DEFINITIONS =================
                        n - number of iterations

             r - scale factor (different scale-factors take                    
            place, r may be vector in this case the lengths 
                    of r and phi must be equal)

        phi - vector of angles in fractal generator (calculated 
        relative to the horizontal axis connecting end points 
                            of generator)
    
                 xb and yb - coordinates of trunk 
                        t - thickness of root in mm

        tr - scale factor for thickness (different scale factors 
             can take place BUT lengths must mach with 
                         lengths of r and phi)

                dim - dimensions of domain in mm 
            pixels - number of pixels along x and y axis        '''

    # Calculating voxel size given pixel and domain size 
    delta_x = dim[0]/pixels[1]
    delta_y = dim[1]/pixels[0]

    # Initializing matrix for domain
    matrix = np.zeros((pixels[1], pixels[0]))

    # Initializnig matrix for scale factors
    mN = len(phi)
    if len(r) == 1:
        rM = np.ones(mN) * r
    elif len(r) == mN:
        rM = np.array(r)
    else:
        print('The sizes of scale vector of lengths and vector of angles in fractal`s generator don`t equal')
        return
    
    if len(tr) == 1:
        tM = int(np.ones(mN) * tr/delta_x)
    elif len(tr) == mN:
         tM = np.array(tr)
    else:
        print('The sizes of scale vector of thickness and vector of angles in fractal`s generator don`t equal')
        return

    # Calculate angle of trunk with vertical axis
    a = np.sqrt((xb[0] /delta_x - xb[1])**2 + (yb[0] - yb[1])**2)
    c = np.sqrt((1 - xb[1])**2 + yb[1]**2)
    alpha = np.arccos((a**2 + 1 - c**2) / (2 * a))
    ksi = alpha if yb[1] >= yb[0] else -alpha

    # Initialize auxiliary vectors
    psi = np.zeros((n, mN**n))
    ralt = np.ones((n, mN**n))
    talt = np.ones((n, mN**n))

    # Assembling matrices through loops
    for i in range(1,n+1):
        z = 1
        for j in range(1, mN**(i-1)+1):
            for k in range(1,mN+1):
                for m in range(1, mN**(n-i)+1):
                    psi[i-1][z-1] = phi[k-1]
                    ralt[i-1][z-1] = rM[k-1]
                    talt[i-1][z-1] = tM[k-1]

                    z += 1

    # Initialize vectors for angles and branch lengths
    theta = np.zeros((n, mN**n))
    rD = np.ones((n, mN**n))
    tD = np.ones((n, mN**n))

    for i in range(1,mN**n+1):
        for j in range(1,n+1):
            for k in range(1,j+1):
                theta[j-1][i-1] += psi[k-1][i-1]
                rD[j-1][i-1] *= ralt[k-1][i-1]
                tD[j-1][i-1] *= talt[k-1][i-1]

    theta += ksi

    # Initialize matrix for coordinates
    A = np.ones((n + 1, mN**n, 2))

    # Initial coordinates
    A[0, :, 0] = np.ones(mN**n) * xb[1]/(delta_x) + pixels[0]//2
    A[0, :, 1] = np.ones(mN**n) * yb[1]/(delta_y)

    # Calculate subsequent coordinates
    for j in range(1,mN**n+1):
        for i in range(1,n+1):
            A[i][j-1][0] = A[i-1][j-1][0] + (a * rD[i-1][j-1] * np.cos(theta[i-1][j-1]))/(delta_x)
            A[i][j-1][1] =  A[i-1][j-1][1] + (a * rD[i-1][j-1] * np.sin(theta[i-1][j-1]))/(delta_y) 

    # Visualizing lines using Bresenhams line algorithm
    for i in range(1, mN**n+1, mN):
        z = 1
        for k in range(1,mN):
            for j in range(z, n ):
                matrix = bresenham2D(A[j-1, i-1, 0], A[j, i-1, 0], A[j-1, i-1, 1], A[j, i-1, 1], int(t*tD[j-1, i-1]*np.cos(theta[j-1][i-1])/delta_x), matrix)
            z += 1

    # Visualizing end branches
    for i in range(mN**n):
        matrix = bresenham2D(A[n-1, i-1, 0], A[n, i-1, 0], A[n-1, i-1, 1], A[n , i-1, 1], int(t*tD[n-1, i-1]*np.cos(theta[n-1][i-1])/delta_x), matrix)

    # Visualizing the trunk
    matrix = bresenham2D(xb[0]/delta_x + pixels[0]//2, xb[1]/delta_x + pixels[0]//2, yb[0]/delta_y, yb[1]/delta_y, int(t*np.sin(alpha)/delta_x), matrix)
    return matrix

matrix = FraktalT(3, [0.8], [np.pi/4, -np.pi/3], [0,0], [0,10], 1, [0.4, 0.8], [50,50], [256,256])
plt.imshow(matrix, origin = 'lower')
plt.show()

