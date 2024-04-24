import numpy as np
import matplotlib.pyplot as plt

# read root from configs directory
def read_txt_to_matrix(file_path):
    matrix = []
    with open(file_path, 'r') as file:
        for line in file:
            element = int(line.strip())  
            matrix.append(element) 
    return matrix

# read domain pixel from configs directory
def read_domain_pixels(file_path):
    with open(file_path, 'r') as config_file:
        for line in config_file:
            if line.startswith("DOMAIN_PIXELS"):
                domain_pixels = line.strip().split('=')[1].strip()[1:-1].split(',')
                domain_pixels = [int(pixel.strip()) for pixel in domain_pixels]
                return domain_pixels
          
# Example usage:
file_path = r'C:\Users\amirt\Desktop\RA\Fractal\configs\ root(1).txt'  
config_file_path = r'C:\Users\amirt\Desktop\RA\Fractal\configs\ config(1).txt'

matrix = read_txt_to_matrix(file_path)
pixels = read_domain_pixels(config_file_path)

# reshaping matrix
reshaped_matrix = np.reshape(matrix, pixels)

# plotting results
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.invert_xaxis()
ax.voxels(reshaped_matrix, facecolor = 'brown', edgecolor = 
              'black', alpha = 0.9, linewidth = 0.5, shade=None)
plt.show() 
