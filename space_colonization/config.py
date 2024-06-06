# config.py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import h5py 
from scipy.signal import convolve
from skan import Skeleton, summarize

# Geometry Options
NO_OF_ITERATIONS = 50                 # NO. OF ITERATIONS THAT SCA RUNS
ROOT_START = [0, 0, 0]               # COORDINATES FOR ORIGIN OF ROOT
ROOT_END = [0, 0, 30]                #  COORDINATES FOR END OF ROOT
ROOT_THICKNESS = 6                   # MAIN THICKNESS OF ROOT
MIN_THICKNESS_BRANCH = 1             # MINIMUM THICKNESS OF ROOT
MAX_THICKNESS_BRANCH = 3             # MAXIMUM THICKNESS OF ROOT
DOMAIN_DIMENSIONS = [20, 20, 60]    # REAL DIMENSIONS OF DOMAIN
DOMAIN_PIXELS = [256, 256, 768]      # NUMBER OF PIXELS ALONG EACH AXIS
BUFFER_LAYER = 0                     # SYMMETRICAL BUFFER LAYER ALONG Z-AXIS

# Space Colonization Options
RADIUS_OF_INFLUENCE = 100                  
KILL_DISTANCE = 1
D = 2
GRAV_ALPHA = 0.4                    # HIGHER VALUE FORCES THE BRANCHES TO GO MORE DOWNARDS
HORIZONTAL_FORCING = 0.4            # HIGHER VALUES FORCES THE BRANCHES TO GO OUT MORE LATERALLY
K_GRAV = 0                          # EXPONENTIAL DECAY CONSTANT FOR THE GRAV_ALPHA OPTION WITH ROOT DEPTH
K_HORIZONTAL = 0                    # EXPONENTIAL DECAY CONSTANT FOR THE HORIZONTAL_FORCING OPTION WITH ROOT DEPTH
         
# Distribution and Root Type
CROWN_TYPE = r'CY'
TAP_ROOT_STYLE = True
no_of_points = 40

#### Cylindrical Settings ('CY') #####
z_min = 0*ROOT_END[2]
z_max = ROOT_END[2]
r_inner = 5
r_outer = 10

#### Ellipsoid Settings ('E') #####
mean = [ROOT_END[0], ROOT_END[1], DOMAIN_DIMENSIONS[2]//2]  
cov = [[DOMAIN_DIMENSIONS[0]**2, 0, 0], 
         [0, DOMAIN_DIMENSIONS[1]**2, 0], 
         [0,0,DOMAIN_DIMENSIONS[2]**2]]
t_x, t_y, t_z = 20, 20, 1
r = DOMAIN_DIMENSIONS[0]

#### Cuboid settings ('C') #####
x_min,x_max = -8, 8
y_min,y_max = -8,8
z_min, z_max = 0, DOMAIN_DIMENSIONS[2]//2

# Saving Options
SAVE_LOCATION = r'C:\Users\...'
FILE_NAME = r'test'
CONFIG_FILE_NAME = r'test_config'
SAVE_AS_TXT = False
SAVE_AS_H5 = True
SAVE_CONFIG_TXT = True

# Other Options
VOLUME_THRESHOLD = 1668356                  # SPECIFY MAXIMUM THRESHOLD OF VOLUME PIXELS OF ROOT
SHOWCASE_RESULT = False
SHOW_ANGULAR_DISTRIBUTION = False
SHOWCASE_BRANCHING_SUMMARY = False
COORDS_OUT_OF_BOUNDS = True
SEED = 581263418
