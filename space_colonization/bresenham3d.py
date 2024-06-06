import numpy as np
from config import DOMAIN_PIXELS, COORDS_OUT_OF_BOUNDS

def bresenham3D(x1, x2, y1, y2, z1, z2, th, matrix):

    voxel_count = 0

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

    # Case I: dx > dy and dx > dz
    if (dx >= dy and dx >= dz): 
            p1 = dy_2 - dx
            p2 = dz_2- dx
            while (x1 != x2):
                
                if y1 >= DOMAIN_PIXELS[1] or z1 >= DOMAIN_PIXELS[2] or x1 >= DOMAIN_PIXELS[0] or x1 < 0 or y1 < 0 or z1 < 0:
                        if not COORDS_OUT_OF_BOUNDS:
                            raise ValueError("Error. Coordinate out of bounds. Terminated code.")
                        break
                else:
                    if int(th) == 1:
                      if x1 >= DOMAIN_PIXELS[0] or x1 < 0 or y1 >= DOMAIN_PIXELS[1] or y1 < 0 or z1 >= DOMAIN_PIXELS[2] or z1 < 0:
                        break
                      else:
                        matrix[y1][x1][z1] = 1
                        voxel_count += 1
                    else:
                        for y_th in range(0, 2*th+1):
                            max_t_z = int(np.sqrt(th**2-(y_th-th)**2))
                            for z_th in range(0,2*max_t_z+1):
                                if y1 - th + y_th >= DOMAIN_PIXELS[0] or y1 - th + y_th < 0 or z1-max_t_z+z_th >= DOMAIN_PIXELS[2] or z1-max_t_z+z_th < 0:
                                   break
                                else:
                                  matrix[y1-th+y_th][x1][z1-max_t_z+z_th] = 1
                                  voxel_count += 1

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
                        if not COORDS_OUT_OF_BOUNDS:
                            raise ValueError("Error. Coordinate out of bounds. Terminated code.")
                        break
                else:
                    if int(th) == 1:
                      if x1 >= DOMAIN_PIXELS[0] or x1 < 0 or y1 >= DOMAIN_PIXELS[1] or y1 < 0 or z1 >= DOMAIN_PIXELS[2] or z1 < 0:
                        break
                      else:
                        matrix[y1][x1][z1] = 1
                        voxel_count += 1 
                    else:
                        for x_th in range(0, 2*th+1):
                            max_t_z = int(np.sqrt(th**2-(x_th-th)**2))
                            for z_th in range(0,2*max_t_z+1):
                                if x1 - th + x_th >= DOMAIN_PIXELS[0] or x1 - th + x_th < 0 or z1 - max_t_z + z_th >= DOMAIN_PIXELS[2] or z1 - max_t_z + z_th < 0:
                                   break
                                else:
                                  matrix[y1][x1-th+x_th][z1-max_t_z+z_th] = 1
                                  voxel_count +=1
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
                        if not COORDS_OUT_OF_BOUNDS:
                            raise ValueError("Error. Coordinate out of bounds. Terminated code.")
                        break
            else:
                if int(th) == 1:
                    if x1 >= DOMAIN_PIXELS[0] or x1 < 0 or y1 >= DOMAIN_PIXELS[1] or y1 < 0 or z1 >= DOMAIN_PIXELS[2] or z1 < 0:
                       break
                    else:
                      matrix[y1][x1][z1] = 1
                      voxel_count += 1
                else:
                    for x_th in range(0, 2*th+1):
                        max_t_y = int(np.sqrt(th**2-(x_th-th)**2))
                        for y_th in range(0,2*max_t_y+1):
                            if y1-max_t_y + y_th >= DOMAIN_PIXELS[0] or  y1-max_t_y + y_th < 0 or x1-th+x_th >= DOMAIN_PIXELS[0] or x1-th+x_th < 0:
                               break
                            else:
                              matrix[y1-max_t_y + y_th][x1-th+x_th][z1] = 1
                              voxel_count +=1 
            if (p1 >= 0):
                y1 += ys
                p1 -= dz_2
            if (p2 >= 0):
                x1 += xs
                p2 -= dz_2
            p1 += dy_2
            p2 += dx_2
            z1 += zs

    return matrix, voxel_count