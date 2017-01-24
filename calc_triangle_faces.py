# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:21:03 2017

@author: Lukas Gartmair
"""

import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import norm
import matplotlib.mlab as mlab

vdb_verts = np.genfromtxt('sharp_vdb_vertices.txt')
vdb_faces = np.genfromtxt('sharp_vdb_faces.txt',dtype=np.int64)
vdb_faces = vdb_faces - 1

ivas_verts = np.genfromtxt('sharp_ivas_vertices.txt')
ivas_faces = np.genfromtxt('sharp_ivas_faces.txt',dtype=np.int64)
ivas_faces = ivas_faces - 1

p12 = np.zeros(3)
p13 = np.zeros(3)

vdb_triangle_areas =  []
for i,f in enumerate(vdb_faces):
    p12[0] = vdb_verts[f[1]][0] - vdb_verts[f[0]][0]
    p12[1] = vdb_verts[f[1]][1] - vdb_verts[f[0]][1]
    p12[2] = vdb_verts[f[1]][2] - vdb_verts[f[0]][2]
    
    p13[0] = vdb_verts[f[2]][0] - vdb_verts[f[1]][0]
    p13[1] = vdb_verts[f[2]][1] - vdb_verts[f[0]][1]
    p13[2] = vdb_verts[f[2]][2] - vdb_verts[f[0]][2]

    crossproduct = np.cross(p12,p13)
    area_parallelogram  = np.linalg.norm(crossproduct)
    
    triangle_area = 0.5 * area_parallelogram   
    
    vdb_triangle_areas.append(triangle_area)
    
    
ivas_triangle_areas =  []
for i,f in enumerate(ivas_faces):
    p12[0] = ivas_verts[f[1]][0] - ivas_verts[f[0]][0]
    p12[1] = ivas_verts[f[1]][1] - ivas_verts[f[0]][1]
    p12[2] = ivas_verts[f[1]][2] - ivas_verts[f[0]][2]
    
    p13[0] = ivas_verts[f[2]][0] - ivas_verts[f[1]][0]
    p13[1] = ivas_verts[f[2]][1] - ivas_verts[f[0]][1]
    p13[2] = ivas_verts[f[2]][2] - ivas_verts[f[0]][2]

    crossproduct = np.cross(p12,p13)
    area_parallelogram  = np.linalg.norm(crossproduct)
    
    triangle_area = 0.5 * area_parallelogram   
    
    ivas_triangle_areas.append(triangle_area)
    
bins = 50

pl.hist(vdb_triangle_areas, bins=bins, label='DMC', color='gray')
pl.grid()
pl.xlabel('triangle area / $nm^2$')
pl.ylabel('frequency')
pl.ylim(0,65)
pl.xlim(0,5)
pl.legend()
pl.figure()
pl.hist(ivas_triangle_areas, bins=bins,label='MC', color='blue')
pl.grid()
pl.xlim(0,5)
pl.ylim(0,65)
pl.xlabel('triangle area / $nm^2$')
pl.ylabel('frequency')
pl.legend()