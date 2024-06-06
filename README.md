# Binary root generation 

The project combines the well-known space colonization algorithm with Bresenhams line algorithm extended to consider thicknesses of lines in order to produce realistic roots and trees. For the ones that wants to understand the algorithms in detail, I heavily recommend this fantastic article [Space Colonization Algorithm](https://ciphrd.com/2019/09/11/generating-a-3d-growing-tree-using-a-space-colonization-algorithm/) by 'ciphrd' and this wikipedia article [Bresenham's Line Algorithm](https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm).

The roots/trees are generated over a binary 3D grid. The intent of this project is to supply researchers with a tool to generate binary voxel representations of roots/trees. These representation can be used as inputs for simulations in computational fluid dynamics and specifically those performed using the lattice Boltzmann method (LBM). 

## Table of contents
- [Installation](#installation)
- [Usage and tips](#usage-and-tips)
- [Examples](#examples)
- [License](#license)

## Installation
In order to run the code, the user will need to supply themselves with the following packages:
- numpy, scipy, matplotlib, h5py, os, skan

## Usage and tips
First configure the parameters in ```config.py```. In the config.py file, one specifies settings such as the root thickness, the number of points sampled in the space colonization algorithm, the domain size in pixels, and much more. There's great variability that can be extracted from this code, so I suggest the user plays around with the parameters themselves.

The ```utils.py``` come with three defined crown types, these are: ellipsoidal, cylindrical and cuboid. The crown types are used for sampling points that are later used in the space colonization algorithm. By choosing different crown types, one can tailor the shape of the root/tree. The user can of course implement new crown types if they are interested. 

The user can also specify whether they want save it as an H5 which has becoming increasingly popular in recent times. For larger domain sizes that are hard to visualize using matplotlib's voxel plot, I suggest the user to install ```VisIt```, see [VisIt installation Guide](https://visit-sphinx-github-user-manual.readthedocs.io/en/v3.2.0/gui_manual/Intro/Installing_VisIt.html).

The second and final step just involves running the main simulation scipt ```main.py```. If wished, one might output the configuration txt file which saves the information from the  ```config.py``` so that one can recreate the exact root/tree in the future. It also contains information such as the number of surface and volume voxels, and surface-to-volume ratio, and other parameters that might be of interest.

## Examples
This chapter shows some examples of roots that have been generated using the provided code.

<p align="center">
   <img src="https://github.com/tresnjo/binary-root-generation/assets/121384892/5c8a2321-b367-4d09-9ebe-58c325fa14e7" alt="Image 1" style="max-width:100%;">
</p>
<p align="center">
        <em>A flower like structure.</em>
</p>
<p align="center">
   <img src="https://github.com/tresnjo/binary-root-generation/assets/121384892/27a07d18-2b71-4d29-87f3-8e39da1747fa" alt="Image 1" style="max-width:100%;">
</p>
<p align="center">
        <em>A fibrous root system reminiscent of grass roots. </em>
</p>

<p align="center">
   <img src="https://github.com/tresnjo/binary-root-generation/assets/121384892/39b86eb7-211d-4960-804a-300e0d03ffa0" alt="Image 1" style="max-width:100%;">
</p>
<p align="center">
        <em>A taproot structure with thick branches.</em>
</p>

<p align="center">
   <img src="https://github.com/tresnjo/binary-root-generation/assets/121384892/d65f5727-88e4-40ee-8a20-49b922ae2b42" alt="Image 1" style="max-width:100%;">
</p>
<p align="center">
        <em>A taproot structure with high gravity forcing of branches.</em>
</p>


## License
