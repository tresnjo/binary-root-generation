from utils import *

if __name__ == "__main__":
    if CROWN_TYPE == 'C':
        run_cuboid() 
    elif CROWN_TYPE == 'E':
        run_ellipsoid()
    elif CROWN_TYPE == 'CY':
        run_cylindrical()
    
