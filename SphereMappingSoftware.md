# Introduction #

Mesh parameterization is central to a broad spectrum of applications.
In our work "Parallel Computation of Spherical Parameterizations for Mesh Analysis" in Computer & Graphics Journal, Volume 35(3), we present a novel approach to spherical mesh parameterization based on an iterative quadratic solver that is efficiently parallelizable on modern massively parallel architectures like GPUs.
Furthermore, we demonstrate the applicability of our approach to real-time feature detection, mesh decomposition and similarity-based 3D object retrieval.

# Software #
![http://archmind.googlecode.com/svn/wiki/images/parameterization_tool.jpg](http://archmind.googlecode.com/svn/wiki/images/parameterization_tool.jpg)
![http://archmind.googlecode.com/svn/wiki/images/parameterization_sphere.jpg](http://archmind.googlecode.com/svn/wiki/images/parameterization_sphere.jpg)

A parallel implementation of the aforementioned algorithm using the archmind geometric kernel can be found in the examples:

[source code](http://code.google.com/p/archmind/source/browse/trunk/examples/sphere_mapping/)

[executable (64bit)](http://archmind.googlecode.com/files/sphere_map_win64.zip)

# Usage #

Usage: sphere\_map --help

Example: sphere\_map --source=Models\suzanne.obj --target=suzanne\_param.obj

Notes:
  1. Please use only closed genus-0 meshes that are centered at the origin
  1. Use only triangulated meshes
  1. Use dense meshes (with at least 100 vertices)