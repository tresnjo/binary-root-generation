
import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D

### Configuration Options ###

ROOT_START = [0,0,0]                                                # Start of root
ROOT_END = [0,0,200]                                                 # End of root
ROOT_THICKNESS = 2
MIN_THICKNESS_BRANCH = 1                                           # Minimum radius of the branch  
DOMAIN_DIMENSIONS = [500,500,500]                                 # Dimensions of the domain 
DOMAIN_PIXELS = [60, 60, 60]                                     # Number of pixels in binary domain            

SEED = 418                                                          # Seed for random fractal generation

METHOD = 'SPACE'                                                  # FRACTAL or SPACE                   
SHOWCASE_RESULT = True                                              # Showcase result
COORDS_OUT_OF_BOUNDS = True                                           # Should coords out of bounds be allowed (T : yes, F: no)
SAVE_MATRIX = False                                                   # Save matrix as HDF5 file
SAVE_CONFIG_TXT = False                                               # Save configuration txt file
SAVE_LOCATION = r'C:\Users\amirt\Desktop\RA\Fractal'                # Saving location
FILE_NAME = r'matrix.h5'                                               # File name

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
    
    th = max(th, MIN_THICKNESS_BRANCH)

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
                        matrix[y1][x1][z1] = 1
                    else:
                        for y_th in range(0, 2*th+1):
                            max_t_z = int(np.sqrt(th**2-(y_th-th)**2))
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
                        if not COORDS_OUT_OF_BOUNDS:
                            raise ValueError("Error. Coordinate out of bounds. Terminated code.")
                        break
                else:
                    if int(th) == 1:
                        matrix[y1][x1][z1] = 1
                    else:
                        for x_th in range(0, 2*th+1):
                            max_t_z = int(np.sqrt(th**2-(x_th-th)**2))
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
                        if not COORDS_OUT_OF_BOUNDS:
                            raise ValueError("Error. Coordinate out of bounds. Terminated code.")
                        break
            else:
                if int(th) == 1:
                    matrix[y1][x1][z1] = 1
                else:
                    for x_th in range(0, 2*th+1):
                        max_t_y = int(np.sqrt(th**2-(x_th-th)**2))
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

    self.d_i = radius_of_influence
    self.d_k = kill_distance
    self.D = D

    self.iter_num = 0

    # attraction points
    x, y, z = crown_attraction_points

    attraction_pts = []
    for i,j,k in list(zip(x, y, z)):
      attraction_pts.append(Attraction_point(i, j,k))

    # nodes
    self.nodes = []
    root = Tree_node(ROOT_END[0], ROOT_END[1], ROOT_END[2])
    self.nodes.append(root)

    # closest node to each attraction pt
    self.closest_node = {attr_pt: None for attr_pt in attraction_pts}
    self.closest_dist = {attr_pt: np.inf for attr_pt in attraction_pts}

    self._update_closest_node(self.nodes[0])

    # branches
    self.branches = []
    self.branching_tree = Tree(root)
    self.brach_min_width = MIN_THICKNESS_BRANCH
    self.branch_width = {}


  def _update_closest_node(self, node):
    kill_candidates = []

    # internal method to update self.closest_node and self.closest_dist
    for attr_pt in self.closest_node:
      old_smallest = self.closest_dist[attr_pt]
      dist = np.linalg.norm(attr_pt.pos - node.pos)

      if dist < self.d_k:
        # attr_pt to be killed
        kill_candidates.append(attr_pt)
        continue

      if dist < self.d_i and dist < old_smallest:
        self.closest_node[attr_pt] = node
        self.closest_dist[attr_pt] = dist

    # kill attraction points with nodes too close to them
    for attr_pt in kill_candidates:
      del self.closest_node[attr_pt]
      del self.closest_dist[attr_pt]
 

  def _report(self):
    print(f'\tREPORT FOR ITERATION {self.iter_num}')
    print('Number of nodes:', len(self.nodes))
    print('Number of attraction points:', len(self.closest_node))
    print()

    data = np.stack([node.pos for node in self.nodes])

    x = []
    y = []
    z = []

    for node in self.nodes:
      x.append(node.x)
      y.append(node.y)
      z.append(node.z)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.scatter(x, y, z, c='slategray')
    plt.show()


  def branch_thinkness(self, node):
    if node in self.branch_width:
      return self.branch_width[node]

    if self.branching_tree.is_leaf(node):
      self.branch_width[node] = self.brach_min_width
      return self.brach_min_width
    
    if self.branching_tree.num_children(node) == 1:
      w = self.branch_thinkness(self.branching_tree.get_children(node)[0])
      self.branch_width[node] = w
      return w

    w = 0
    for child in self.branching_tree.get_children(node):
      w += np.square(self.branch_thinkness(child))
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


  def render_results(self):
    
    matrix = np.zeros((DOMAIN_PIXELS[0], DOMAIN_PIXELS[1], DOMAIN_PIXELS[2]))

    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    delta_x = DOMAIN_DIMENSIONS[0]/DOMAIN_PIXELS[0]
    delta_y = DOMAIN_DIMENSIONS[1]/DOMAIN_PIXELS[1]
    delta_z = DOMAIN_DIMENSIONS[2]/DOMAIN_PIXELS[2]

    x = np.array([ROOT_START[0], ROOT_END[0]])/delta_x + DOMAIN_PIXELS[0]//2
    y = np.array([ROOT_START[1], ROOT_END[1]])/delta_y+ DOMAIN_PIXELS[1]//2
    z = np.array([ROOT_START[2], ROOT_END[2]])/delta_z 

      #ax.plot(x, y, z, c='forestgreen', linewidth=lw)
    matrix = bresenham3D(x[0], x[1], y[0], y[1], z[0], z[1], 3, matrix)

    for branch in self.branches:
      start, end, node = branch

      lw = self.branch_thinkness(node)

      x = np.array([start[0], end[0]])/delta_x + DOMAIN_PIXELS[0]//2
      y = np.array([start[1], end[1]])/delta_y+ DOMAIN_PIXELS[1]//2
      z = np.array([start[2], end[2]])/delta_z 

      #ax.plot(x, y, z, c='forestgreen', linewidth=lw)
      matrix = bresenham3D(x[0], x[1], y[0], y[1], z[0], z[1], int(round(lw/4)), matrix)

    ax.invert_xaxis()
    ax.set_axis_off()
    ax.voxels(matrix, facecolor = 'brown', edgecolor = 
            'black', alpha = 0.9, linewidth = 0.5, shade=None)
    plt.show()


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
          n += (attr_pt.pos - node.pos) / np.linalg.norm(attr_pt.pos - node.pos)
        n = n / np.linalg.norm(n)

        new_pos = node.pos + n * self.D
        new_node = Tree_node(new_pos[0], new_pos[1], new_pos[2])
        self._update_closest_node(new_node)

        branch = (node.pos, new_pos, new_node)
        self.branches.append(branch)
        self.branching_tree.add_child(node, new_node)

        meta_nodes.append(new_node)
    # add newly added nodes
    self.nodes.extend(meta_nodes)

def run_experiment_ellispe_crown_1():
  
  np.random.seed(SEED)
  
  fig = plt.figure()
  ax = fig.add_subplot(projection='3d')

  mean = [ROOT_END[0], ROOT_END[1], ROOT_END[2]]

  # när du kmr hem fixa kovariansmatrisen så du får en bra spridning och sedan även fixa t saken
  # efter det kan du nog få ut relativt fina resultat
  # justera på parametrarna vidare och se om du kan komma på ngn generell formel
  cov = [[DOMAIN_DIMENSIONS[0]**2, 0, 0], 
         [0, DOMAIN_DIMENSIONS[1]**2, 0], 
         [0,0,DOMAIN_DIMENSIONS[2]**2]]

  x, y, z= np.random.multivariate_normal(mean, cov, 1000).T

  t = 3*np.square(x-ROOT_END[0])+ 3*np.square(y-ROOT_END[1]) + 1.5*np.square(z-ROOT_END[2]) <= DOMAIN_DIMENSIONS[0]**2

  x_crown = x[t]
  y_crown = y[t]
  z_crown = z[t]

  ax.plot(x_crown, y_crown, z_crown, 'o')
  plt.show()

  sim = Simulation(crown_attraction_points=(x_crown, y_crown, z_crown), radius_of_influence = 200, kill_distance = 1, D = 10)
  sim.run(50)

  del sim

run_experiment_ellispe_crown_1() 
     
