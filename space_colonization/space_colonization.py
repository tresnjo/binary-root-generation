from tree import *
from bresenham3d import bresenham3D
from config import *
from utils import *

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
    self.filter_isolated_islands = np.array([
                              [[1, 1, 1],    
                               [1, 1, 1],
                               [1, 1, 1]],

                              [[1, 1, 1],
                               [1 ,0, 1],
                               [1, 1, 1]],

                              [[1, 1, 1],
                               [1, 1, 1],
                               [1, 1, 1]]
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
    
  def count_isolated_islands(self, matrix):
     mask = convolve(matrix, self.filter_isolated_islands, mode = 'same')
     pixels = (matrix == 1) & (mask == 0)
     self.isolated_islands = np.sum(pixels)
  
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

  # returns a summary of the branches type, lengths etc
  def branching_summary(self, root_matrix, delta_x, delta_y, delta_z):
    ''' 
    Notes on histogram:
        0: endpoint to endpoint
        1: junction to endpoint
        2: junction to junction
        3: isolated cycle      
    '''
    branch_data = summarize(Skeleton(root_matrix, spacing = [delta_x, delta_y, delta_z]))
    if SHOWCASE_BRANCHING_SUMMARY:
      branch_data.hist(column='branch-distance', by='branch-type', bins=50)
      plt.show()
    return branch_data
     
  def render_results(self):
    
    # initializing matrix from specified pixel dimensions
    matrix = np.zeros((DOMAIN_PIXELS[0], DOMAIN_PIXELS[1], DOMAIN_PIXELS[2])).astype(float)
    pseudo_matrix = np.zeros((DOMAIN_PIXELS[0], DOMAIN_PIXELS[1], DOMAIN_PIXELS[2])).astype(float)
  
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
      max_dev_x = 15
      max_dev_y = 15
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
        matrix, new_pixels = bresenham3D(x_init[0], x_init[1], y_init[0], y_init[1], z_init[0], 0, int(ROOT_THICKNESS), matrix)
        #matrix, new_pixels = bresenham3D(x_init[0], x_init[1], y_init[0], y_init[1], z_init[0], z_init[1], int(ROOT_THICKNESS), matrix)
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
        self.acc_pixels += new_pixels
      else:
        print("\nVolume threshold exceeded.")
        break

    # calculating number of surface voxels
    self.count_surface_area(matrix)
    self.count_isolated_islands(matrix)
    branching_data = self.branching_summary(pseudo_matrix, delta_x, delta_y, delta_z)

    # adding buffer layer
    buffer_matrix = np.zeros((DOMAIN_PIXELS[0], DOMAIN_PIXELS[1], BUFFER_LAYER))
    matrix = np.concatenate((matrix, buffer_matrix), axis = 2)

    # inverting matrix
    matrix = matrix[:,:,::-1]
    pseudo_matrix = pseudo_matrix[:,:,::-1]

    # adding second buffer layer
    matrix = np.concatenate((matrix, buffer_matrix), axis = 2)

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
            config_file.write("K_GRAVITY = {}\n".format(K_GRAV))
            config_file.write("MIN_THICKNESS_BRANCH = {}\n".format(MIN_THICKNESS_BRANCH))
            config_file.write("MAX_THICKNESS_BRANCH = {}\n".format(MAX_THICKNESS_BRANCH))
            config_file.write("DOMAIN_DIMENSIONS = {}\n".format(DOMAIN_DIMENSIONS))
            config_file.write("DOMAIN_PIXELS = {}\n".format(DOMAIN_PIXELS))
            config_file.write("BUFFER_LAYER = {}\n".format(BUFFER_LAYER))
            config_file.write("SEED = {}\n".format(SEED))
            config_file.write("RADIUS_OF_INFLUENCE = {}\n".format(RADIUS_OF_INFLUENCE))
            config_file.write("KILL_DISTANCE = {}\n".format(KILL_DISTANCE))
            config_file.write("D = {}\n".format(D))
            config_file.write("TAP_ROOT_STYLE = {}\n".format(TAP_ROOT_STYLE))
            config_file.write("HORIZONTAL_FORCING = {}\n".format(HORIZONTAL_FORCING))
            config_file.write("K_HORIZONTAL = {}\n".format(K_HORIZONTAL))
            config_file.write("VOLUME_THRESHOLD = {}\n".format(VOLUME_THRESHOLD))

            config_file.write("\n#### ROOT RESULTS #### \n")
            config_file.write("VOL_PIXELS = {}\n".format(self.acc_pixels))
            config_file.write("SURF_VOL_PIXELS = {}\n".format(self.surf_pixels))
            config_file.write("POROSITY = {}\n".format(self.acc_pixels/(DOMAIN_PIXELS[0]*DOMAIN_PIXELS[1]*DOMAIN_PIXELS[2])))
            config_file.write("SURF_TO_VOL_RATIO = {}\n".format(self.surf_pixels/self.acc_pixels))
            config_file.write("NO_OF_BRANCHES = {}\n".format(branching_data.shape[0]))
            config_file.write("NO_OF_ISOLATED_PIXELS = {}\n".format(self.isolated_islands))

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
          hf.create_dataset("branching_data", data=branching_data)
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
          delta_x_sgn = np.sign(attr_pt.pos[0])
          delta_y_sgn = np.sign(attr_pt.pos[1])
          n += n_c + GRAV_ALPHA*np.array([0,0,1])*np.exp(K_GRAV*attr_pt.pos[2]/ROOT_END[2]) + HORIZONTAL_FORCING*np.array([delta_x_sgn*np.random.rand(), delta_y_sgn*np.random.rand(), 0])*np.exp(-K_HORIZONTAL*attr_pt.pos[2]/ROOT_END[2])

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
