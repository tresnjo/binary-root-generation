from config import *
from simulation import *

def run_cuboid():
  np.random.seed(SEED)

  global x_min, x_max, y_min, y_max, z_min, z_max, no_of_points

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

  x, y, z= np.random.multivariate_normal(mean, cov, no_of_points).T

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