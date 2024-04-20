
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import h5py 
from scipy.signal import convolve
from skan import Skeleton, summarize

########################################################
#########       CONFIGURATION OPTIONS       ############
########################################################

''' GEOMETRY OPTIONS'''
NO_OF_ITERATIONS = 50                                             # No iterations
PARTICLE_SIZE = 5                                                # Specify in relation to pixels
ROOT_START = [0,0,0]                                                # Start of main root
ROOT_END = [0,0,40]                                                 # End of main root in terms
ROOT_THICKNESS = 6*PARTICLE_SIZE                                                  # Speficication of main root radius
MIN_THICKNESS_BRANCH = 1.5*PARTICLE_SIZE                                          # Minimum radius of the branch  
MAX_THICKNESS_BRANCH = 4*PARTICLE_SIZE                                          # Maximum radius of the branch
DOMAIN_DIMENSIONS = [20,20,80]                                 # Dimensions of the domain 
DOMAIN_PIXELS = [256, 256, 1024]                                     # Number of pixels in binary domain            

''' SPACE COLONIZATION OPTIONS'''
RADIUS_OF_INFLUENCE = 100                                            # Radius of incluence for space colonization algorithm
KILL_DISTANCE = 1                                                # Kill distance for space colonization algorithm
D = 2                                                           # Jump distance D
GRAV_ALPHA = 1                                                    # Gravitropism 
HORIZONTAL_FORCING = 0                                            # Horizontal forcing

''' DISTRIBUTION AND ROOT TYPE'''
CROWN_TYPE = r'CY'                                                # CUBOID (C) or ELLIPSOIDE (E) or CYLINDRICAL (CY)
TAP_ROOT_STYLE = True                                            # To generate tap root of main root

""" SAVING OPTIONS """
SAVE_LOCATION = r'C:\Users\amirt\Desktop\RA\Fractal\Configs\taproots\ '                # Saving location
FILE_NAME = r'test_surf_vol'                                               # Root txt and h5 file name
CONFIG_FILE_NAME = r'test_surf_vol'                                      # Configuration txt file name
SAVE_AS_TXT = False                                                   # Save matrix as binary txt file
SAVE_AS_H5 = False                                                     # Save matrix, polar and azimuthal as H5 file 
SAVE_CONFIG_TXT = True                                               # Save configuration settings as txt file

''' OTHER OPTIONS'''
VOLUME_THRESHOLD = 1e8                                            # Maximum number of volume voxels 
SHOWCASE_RESULT = True                                              # Showcase result (T : yes, F: no)
SHOW_ANGULAR_DISTRIBUTION = False                                   # Showing angular distribution of root
SHOWCASE_BRANCHING_SUMMARY = True                                   # Showing branch length distribution of root
COORDS_OUT_OF_BOUNDS = True                                           # Should coords out of bounds be allowed (T : yes, F: no)
SEED = 581263418                                                          # Seed for random root generation  

########################################################
#######       BRESENHAMS 3D LINE ALGORITHM       #######
########################################################

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

########################################################
#######       SPACE COLONIZATION ALGORITHM       #######
########################################################

class Tree_node:
  def __init__(self, pos_x, pos_y, pos_z):
    self.x = pos_x
    self.y = pos_y
    self.z = pos_z
    self.pos = np.array([self.x, self.y, self.z])

  def __repr__(self):
    return self.pos.__repr__()


class Attraction_point:
  def __init__(self, x, y, z):
    self.x = x
    self.y = y
    self.z = z
    self.pos = np.array([self.x, self.y, self.z])

  def __repr__(self):
    return self.pos.__repr__()
   

class Tree:
  def __init__(self, root):
    self.root = root
    self.nodes = [self.root]
    self.transition_map = {}

  def add_child(self, parent, child):

    if child in self.nodes:
      raise ValueError

    if parent in self.nodes:
      self.transition_map[child] = parent
      self.nodes.append(child)
    else:
      raise ValueError

  def is_leaf(self, node):
    if node not in self.nodes:
      raise ValueError
    if node not in self.transition_map.values():
      return True
    else:
      return False

  def num_children(self, node):
    if node not in self.nodes:
      raise ValueError

    return list(self.transition_map.values()).count(node)

  def get_children(self, parent):
    if parent not in self.nodes:
      raise ValueError

    return [child for child in self.transition_map if self.transition_map[child] == parent]


  def get_level(self, node):
    if node not in self.nodes:
      raise ValueError

    if node == self.root:
      return 0

    x = self.transition_map[node]
    level = 1
    while x != self.root:
      x = self.transition_map[x]

      if self.num_children(x) > 1:
        level += 1

    return level


class Simulation:
  def __init__(self, crown_attraction_points, radius_of_influence, kill_distance, D):

    # variable initialization
    self.acc_pixels = 0
    self.surf_pixels = 0
    self.no_of_branches = 0
    self.thetas = []
    self.phis = []
    self.filter3d = np.array([
                              [[0, 0, 0],    
                               [0, 1, 0],
                               [0, 0, 0]],

                              [[0, 1, 0],
                               [1 ,0, 1],
                               [0, 1, 0]],

                              [[0, 0, 0],
                               [0, 1, 0],
                               [0, 0, 0]]
                              ]
                           )

    self.d_i = radius_of_influence
    self.d_k = kill_distance
    self.D = D
    self.iter_num = 0

    x, y, z = crown_attraction_points

    attraction_pts = []
    for i,j,k in list(zip(x, y, z)):
      attraction_pts.append(Attraction_point(i, j,k))

    # nodes
    self.nodes = []
    root = Tree_node(0, 0, 0)
    self.nodes.append(root)

    # closest node to each attraction pt
    self.closest_node = {attr_pt: None for attr_pt in attraction_pts}
    self.closest_dist = {attr_pt: np.inf for attr_pt in attraction_pts}

    self._update_closest_node(self.nodes[0])

    # branches
    self.branches = []
    self.branching_tree = Tree(root)
    self.branch_min_width = MIN_THICKNESS_BRANCH
    self.branch_width = {}

  # calculating the surface area using convolution of 3D filter
  def count_surface_area(self, matrix):
      mask = convolve(matrix, self.filter3d, mode = 'same')
      boundaries = np.logical_and(mask > 0, mask < 6)
      self.surf_pixels = np.sum(matrix[boundaries])
  
  # updating the closest nodes
  def _update_closest_node(self, node):
    kill_candidates = []

    for attr_pt in self.closest_node:
      old_smallest = self.closest_dist[attr_pt]
      dist = np.linalg.norm(attr_pt.pos - node.pos)

      if dist < self.d_k:
        kill_candidates.append(attr_pt)
        continue

      if dist < self.d_i and dist < old_smallest:
        self.closest_node[attr_pt] = node
        self.closest_dist[attr_pt] = dist

    for attr_pt in kill_candidates:
      del self.closest_node[attr_pt]
      del self.closest_dist[attr_pt]
 

  # report function 
  def _report(self):
    print(f'\tREPORT FOR ITERATION {self.iter_num}')
    print('Number of nodes:', len(self.nodes))
    print('Number of attraction points:', len(self.closest_node))
    print()

    x = []
    y = []
    z = []

    for node in self.nodes:
      x.append(node.x)
      y.append(node.y)
      z.append(node.z)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_box_aspect((1,DOMAIN_DIMENSIONS[1]//DOMAIN_DIMENSIONS[0],DOMAIN_DIMENSIONS[2]//DOMAIN_DIMENSIONS[0]))
    ax.scatter(x, y, z, c='slategray')
    plt.show()

  # getting the branch thickness 
  def branch_thickness(self, node):
    if node in self.branch_width:
      return self.branch_width[node]

    if self.branching_tree.is_leaf(node):
      self.branch_width[node] = self.branch_min_width
      return self.branch_min_width
    
    if self.branching_tree.num_children(node) == 1:
      w = self.branch_thickness(self.branching_tree.get_children(node)[0])
      self.branch_width[node] = w
      return w

    w = 0
    for child in self.branching_tree.get_children(node):
      w += np.square(self.branch_thickness(child))
    w = np.sqrt(w)

    self.branch_width[node] = w
    return w

  def run(self, num_iteration):

    for i in range(num_iteration):
      self._iter()

      if len(self.closest_node) == 0:
        break
    
    self._report()
    self.render_results()

  def branching_summary(self, root_matrix, delta_x):
    ''' Notes _________________
        0: endpoint to endpoint
        1: junction to endpoint
        2: junction to junction
        3: isolated cycle      
    '''
    branch_data = summarize(Skeleton(root_matrix, spacing = delta_x))
    if SHOWCASE_BRANCHING_SUMMARY:
      branch_data.hist(column='branch-distance', by='branch-type', bins=100)
      plt.show()
    return branch_data
     
  def render_results(self):
    
    # initializing matrix from specified pixel dimensions
    matrix = np.zeros((DOMAIN_PIXELS[0], DOMAIN_PIXELS[1], DOMAIN_PIXELS[2])).astype(int)
    pseudo_matrix = np.zeros((DOMAIN_PIXELS[0], DOMAIN_PIXELS[1], DOMAIN_PIXELS[2])).astype(int)
  

    # conversion from domain dimension to pixel dimension
    delta_x = DOMAIN_DIMENSIONS[0]/DOMAIN_PIXELS[0]
    delta_y = DOMAIN_DIMENSIONS[1]/DOMAIN_PIXELS[1]
    delta_z = DOMAIN_DIMENSIONS[2]/DOMAIN_PIXELS[2]

    # initial position from root to pixel dimension
    x_init = np.array([ROOT_START[0], ROOT_END[0]])/delta_x + DOMAIN_PIXELS[0]//2
    y_init = np.array([ROOT_START[1], ROOT_END[1]])/delta_y+ DOMAIN_PIXELS[1]//2
    z_init = np.array([ROOT_START[2], ROOT_END[2]])/delta_z 

    # loop for calculating maximum thickness used later for rescaling
    max_th = 0
    for branch in self.branches:
       start, end, node = branch
       lw = self.branch_thickness(node)
       if lw > max_th:
          max_th = lw
    
    ### TAPROOT STYLE FOR MAIN ROOT ###
    if TAP_ROOT_STYLE: 
      global max_dev_x, max_dev_y, no_of_division, slope

      # number of divisions along main root alongside max deviations as random perturbations in x and y
      # slope is used for specifying how fast thickness of root should decrease 
      max_dev_x = 20
      max_dev_y = 20
      no_of_division = 15
      slope = 0.08

      z = np.linspace(z_init[0], z_init[1], no_of_division)
      xs = np.random.uniform(-max_dev_x, max_dev_x, no_of_division)
      ys = np.random.uniform(-max_dev_y, max_dev_y, no_of_division)

      for i in range(1,no_of_division):
        if self.acc_pixels < VOLUME_THRESHOLD:
          matrix, new_pixels = bresenham3D(x_init[0]+xs[i-1], x_init[1]+xs[i], y_init[0]+ys[i-1], y_init[1]+ys[i], z[i-1], z[i], int(self.tap_root(ROOT_THICKNESS,slope,z[i-1])), matrix)
          self.acc_pixels += new_pixels
        else:
           print("\nVolume threshold exceeded.")
           break
    else:
      if self.acc_pixels < VOLUME_THRESHOLD:
        matrix, new_pixels = bresenham3D(x_init[0], x_init[1], y_init[0], y_init[1], z_init[0], z_init[1], int(ROOT_THICKNESS), matrix)
        self.acc_pixels += new_pixels
      else:
        print("\nVolume threshold exceeded.")

    for branch in self.branches:
      start, end, node = branch

      lw = self.branch_thickness(node)
      x = np.array([start[0], end[0]])/delta_x + DOMAIN_PIXELS[0]//2
      y = np.array([start[1], end[1]])/delta_y+ DOMAIN_PIXELS[1]//2
      z = np.array([start[2], end[2]])/delta_z 

      if self.acc_pixels < VOLUME_THRESHOLD:
        matrix, new_pixels = bresenham3D(x[0], x[1], y[0], y[1], z[0], z[1], int(MIN_THICKNESS_BRANCH + (MAX_THICKNESS_BRANCH - MIN_THICKNESS_BRANCH)*lw/max_th), matrix)
        pseudo_matrix, pseudo_pixels = bresenham3D(x[0], x[1], y[0], y[1], z[0], z[1], 1, pseudo_matrix)
        self.no_of_branches += 1
        self.acc_pixels += new_pixels
      else:
        print("\nVolume threshold exceeded.")
        break
    
    # inverting matrix for root generation
    matrix = matrix[:,:,::-1]
    pseudo_matrix = pseudo_matrix[:,:,::-1]

    # calculating number of surface voxels
    self.count_surface_area(matrix)
    branching_data = self.branching_summary(pseudo_matrix, delta_x)

    # saving the matrix as a binary txt file
    if SAVE_AS_TXT:
      save_path = self.save_file(f"{SAVE_LOCATION}{FILE_NAME}.txt", FILE_NAME, type = "txt")
      np.savetxt(save_path, matrix.flatten().astype(int), fmt = '%.0f')
      print("\nMatrix txt-file successfully created at {}".format(save_path))

    # writing the configuration file
    if SAVE_CONFIG_TXT:
      config_file_path = self.save_file(f"{SAVE_LOCATION}{CONFIG_FILE_NAME}.txt", CONFIG_FILE_NAME, type = "txt")
      with open(config_file_path, 'w') as config_file:
            
            config_file.write("\n#### CONFIGURATION SETTINGS #### \n")
            config_file.write("ROOT_START = {}\n".format(ROOT_START))
            config_file.write("ROOT_END = {}\n".format(ROOT_END))
            config_file.write("ROOT_THICKNESS = {}\n".format(ROOT_THICKNESS))
            config_file.write("GRAVITROPISM = {}\n".format(GRAV_ALPHA))
            config_file.write("MIN_THICKNESS_BRANCH = {}\n".format(MIN_THICKNESS_BRANCH))
            config_file.write("MAX_THICKNESS_BRANCH = {}\n".format(MAX_THICKNESS_BRANCH))
            config_file.write("DOMAIN_DIMENSIONS = {}\n".format(DOMAIN_DIMENSIONS))
            config_file.write("DOMAIN_PIXELS = {}\n".format(DOMAIN_PIXELS))
            config_file.write("SEED = {}\n".format(SEED))
            config_file.write("RADIUS_OF_INFLUENCE = {}\n".format(RADIUS_OF_INFLUENCE))
            config_file.write("KILL_DISTANCE = {}\n".format(KILL_DISTANCE))
            config_file.write("D = {}\n".format(D))
            config_file.write("TAP_ROOT_STYLE = {}\n".format(TAP_ROOT_STYLE))
            config_file.write("HORIZONTAL_FORCING = {}\n".format(HORIZONTAL_FORCING))
            config_file.write("VOLUME_THRESHOLD = {}\n".format(VOLUME_THRESHOLD))

            config_file.write("\n#### ROOT RESULTS #### \n")
            config_file.write("VOL_PIXELS = {}\n".format(self.acc_pixels))
            config_file.write("SURF_VOL_PIXELS = {}\n".format(self.surf_pixels))
            config_file.write("POROSITY = {}\n".format(self.acc_pixels/(DOMAIN_PIXELS[0]*DOMAIN_PIXELS[1]*DOMAIN_PIXELS[2])))
            config_file.write("SURF_TO_VOL_RATIO = {}\n".format(self.surf_pixels/self.acc_pixels))
            config_file.write("NO_OF_BRANCHES = {}\n".format(self.no_of_branches))

            if TAP_ROOT_STYLE:
              config_file.write("\n#### TAPROOT SETTINGS #### \n")
              config_file.write("MAX_DELTA_X = {}\n".format(max_dev_x))
              config_file.write("MAX_DELTA_Y = {}\n".format(max_dev_y))
              config_file.write("NO_OF_DIVS = {}\n".format(no_of_division))
              config_file.write("SLOPE = {}\n".format(slope))
            
            config_file.write("\n#### DISTRIBUTION SETTINGS #### \n")
            
            if CROWN_TYPE == 'C':
              config_file.write("CROWN_TYPE = {}\n".format(CROWN_TYPE))
              config_file.write("x_min, x_max = {}\n".format([x_min,x_max]))
              config_file.write("y_min, y_max = {}\n".format([y_min,y_max]))
              config_file.write("z_min, z_max = {}\n".format([z_min,z_max]))
              config_file.write("no_of_points = {}\n".format(no_of_points))
            elif CROWN_TYPE == 'E':
              config_file.write("CROWN_TYPE = {}\n".format(CROWN_TYPE))
              config_file.write("mean = {}\n".format(mean))
              config_file.write("cov = {}\n".format(cov))
              config_file.write("t_x, t_y, t_z = {}\n".format([t_x,t_y,t_z]))
              config_file.write("no_of_points = {}\n".format(no_of_points))
              config_file.write("r = {}\n".format(r))
            elif CROWN_TYPE == 'CY':
              config_file.write("CROWN_TYPE = {}\n".format(CROWN_TYPE))
              config_file.write("no_of_points = {}\n".format(no_of_points))
              config_file.write("z_min, z_max = {}\n".format([z_min,z_max]))
              config_file.write("r_inner, r_outer = {}\n".format([r_inner,r_outer]))

      print("\nConfiguration file successfully created at {}".format(config_file))
    
    # saving the matrix as a binary hdf5 file
    if SAVE_AS_H5:
      save_path_h5 = self.save_file(f"{SAVE_LOCATION}{FILE_NAME}.h5", FILE_NAME, type = "h5")
      with h5py.File(save_path_h5, 'w') as hf:
          hf.create_dataset("root",  data=matrix)
          hf.create_dataset("branching_data", data=pseudo_matrix)
          hf.create_dataset("azimuth", data = self.phis)
          hf.create_dataset("polar", data = self.thetas)
          hf['root'].attrs['element_size_mm'] = [delta_x, delta_y, delta_z]
          
      print("\nH5-file successfully created at {}".format(save_path_h5))

    # showcasing the matrix representation of the root using voxels
    if SHOWCASE_RESULT:
      # initializing plot
      fig = plt.figure()
      ax = fig.add_subplot(projection='3d')
      ax.set_box_aspect((1,DOMAIN_DIMENSIONS[1]//DOMAIN_DIMENSIONS[0],DOMAIN_DIMENSIONS[2]//DOMAIN_DIMENSIONS[0]))
      ax.invert_xaxis()
      ax.voxels(matrix, facecolor = 'brown', edgecolor = 
              'black', alpha = 0.9, linewidth = 0.5, shade=None)
      ax.set_box_aspect((1,DOMAIN_DIMENSIONS[1]//DOMAIN_DIMENSIONS[0],DOMAIN_DIMENSIONS[2]//DOMAIN_DIMENSIONS[0]))
      plt.show()
    
    # showcasing the angular distributions of root
    if SHOW_ANGULAR_DISTRIBUTION:
      plt.style.use("seaborn-v0_8-paper")
      fig, axs = plt.subplots(1, 2, figsize=(12, 6), subplot_kw=dict(projection='polar'))
      fig.suptitle("Branching angle distribution", fontsize = 16)
      axs[0].hist(self.phis, bins=30, color='skyblue', alpha=0.3, density = True)
      axs[0].set_title("Distribution of azimuthal angle " + r"$\psi$")
      axs[1].hist(self.thetas, bins=30, color='skyblue', alpha=0.3, density = True)
      axs[1].set_title("Distribution of polar angle " + r"$\theta$")
      plt.show()

    return  
  
  # taproot linear function for thickness as function of depth
  def tap_root(self, m, k, z):
     return m - k * z
  
  # iterating over saving paths
  def save_file(self,save_path, filename, type):
    if os.path.exists(save_path):
        i = 1
        while True:
          new_filename = f"{SAVE_LOCATION}{filename}({i}).{type}"
          if not os.path.exists(new_filename):
              save_path = new_filename
              break
          i += 1
    return save_path
  
  def _iter(self):

    self.iter_num += 1

    meta_nodes = []
    for node in self.nodes:
      # find set of attraction pts affecting node
      S_v = {attr_pt for attr_pt in self.closest_node if self.closest_node[attr_pt] == node}

      # if set is not empty, add new node
      if len(S_v) != 0:
        # find new node pos
        n = np.array([0, 0, 0], dtype=float)
        for attr_pt in S_v:
          n_c = (attr_pt.pos - node.pos) / np.linalg.norm(attr_pt.pos - node.pos)
          n += n_c + GRAV_ALPHA*np.array([0,0,1]) + HORIZONTAL_FORCING*np.array([np.random.rand(), np.random.rand(), 0])

        n = n / np.linalg.norm(n)

        self.phis.append(np.arctan2(n[1],n[0])*180/np.pi)
        self.thetas.append(np.arccos(n[2]/np.linalg.norm(n))*180/np.pi)

        new_pos = node.pos + n * self.D
        new_node = Tree_node(new_pos[0], new_pos[1], new_pos[2])
        self._update_closest_node(new_node)

        branch = (node.pos, new_pos, new_node)
        self.branches.append(branch)
        self.branching_tree.add_child(node, new_node)

        meta_nodes.append(new_node)
    # add newly added nodes
    self.nodes.extend(meta_nodes)


########################################################
#######       VARIOUS SPACE DISTRIBUTIONS       ########
########################################################
def run_cuboid():
  np.random.seed(SEED)

  global x_min, x_max, y_min, y_max, z_min, z_max, no_of_points
  no_of_points = 400
  x_min = -8
  x_max = 8
  y_min = -8
  y_max = 8
  z_min = 0
  z_max = DOMAIN_DIMENSIONS[2]//2

  x_coords = np.random.uniform(x_min, x_max, no_of_points)
  y_coords = np.random.uniform(y_min, y_max, no_of_points)
  z_coords = np.random.uniform(z_min, z_max, no_of_points)

  sim = Simulation(crown_attraction_points=(x_coords, y_coords, z_coords), radius_of_influence = RADIUS_OF_INFLUENCE, kill_distance = KILL_DISTANCE, D = D)
  sim.run(NO_OF_ITERATIONS)

def run_ellipsoid():
  
  np.random.seed(SEED)
  
  fig = plt.figure()
  ax = fig.add_subplot(projection='3d')

  global mean, cov, t_x, t_y, t_z, no_of_points, r

  no_of_points = 2000

  mean = [ROOT_END[0], ROOT_END[1], DOMAIN_DIMENSIONS[2]//2]
  
  cov = [[DOMAIN_DIMENSIONS[0]**2, 0, 0], 
         [0, DOMAIN_DIMENSIONS[1]**2, 0], 
         [0,0,DOMAIN_DIMENSIONS[2]**2]]

  x, y, z= np.random.multivariate_normal(mean, cov, no_of_points).T

  t_x = 20
  t_y = 20
  t_z = 1
  r = DOMAIN_DIMENSIONS[0]

  t = t_x*np.square(x-mean[0])+ t_y*np.square(y-mean[1]) + t_z * np.square(z-(mean[2])) <= r**2

  x_crown = x[t]
  y_crown = y[t]
  z_crown = z[t]

  ax.plot(x_crown, y_crown, z_crown, 'o')
  ax.set_box_aspect((1,DOMAIN_DIMENSIONS[1]//DOMAIN_DIMENSIONS[0],DOMAIN_DIMENSIONS[2]//DOMAIN_DIMENSIONS[0]))
  plt.show()

  sim = Simulation(crown_attraction_points=(x_crown, y_crown, z_crown), radius_of_influence = RADIUS_OF_INFLUENCE, kill_distance= KILL_DISTANCE, D = D)
  sim.run(NO_OF_ITERATIONS)

  del sim
   
def run_cylindrical():

  np.random.seed(SEED)

  fig = plt.figure()
  ax = fig.add_subplot(projection='3d')
    
  global no_of_points, z_min, z_max, r_inner, r_outer

  no_of_points = 150
  z_min = 0*ROOT_END[2]
  z_max = ROOT_END[2]
  r_inner = 4#0*ROOT_THICKNESS/DOMAIN_PIXELS[0] * DOMAIN_DIMENSIONS[0]
  r_outer = 12 #3*ROOT_THICKNESS/DOMAIN_PIXELS[0] * DOMAIN_DIMENSIONS[0]
  r = np.random.uniform(r_inner, r_outer, no_of_points)
  theta = np.random.uniform(0, 2*np.pi, no_of_points)
  z = np.random.uniform(z_min, z_max, no_of_points)

  x =  r*np.cos(theta)
  y = r*np.sin(theta)

  ax.plot(x, y, z, 'o')
  ax.set_box_aspect((1,DOMAIN_DIMENSIONS[1]//DOMAIN_DIMENSIONS[0],DOMAIN_DIMENSIONS[2]//DOMAIN_DIMENSIONS[0]))
  plt.show()

  sim = Simulation(crown_attraction_points=(x, y, z), radius_of_influence = RADIUS_OF_INFLUENCE, kill_distance = KILL_DISTANCE, D = D)
  sim.run(NO_OF_ITERATIONS)
  

if CROWN_TYPE == 'C':
  run_cuboid() 
elif CROWN_TYPE == 'E':
  run_ellipsoid()
elif CROWN_TYPE == 'CY':
  run_cylindrical()

