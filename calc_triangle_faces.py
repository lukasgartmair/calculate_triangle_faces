# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:21:03 2017

@author: Lukas Gartmair
"""

import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import norm
import matplotlib.mlab as mlab

vdb_verts = np.genfromtxt('erboco_vdb_vertices.txt')
vdb_faces = np.genfromtxt('erboco_vdb_faces.txt',dtype=np.int64)
vdb_faces = vdb_faces - 1

ivas_verts = np.genfromtxt('erboco_ivas_vertices.txt')
ivas_faces = np.genfromtxt('erboco_ivas_faces.txt',dtype=np.int64)
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

fig, ax = pl.subplots()

ax.hist(vdb_triangle_areas, bins=bins, label='DMC', color='gray')
pl.grid()
pl.tick_params(axis='both', which='major', labelsize=19)
pl.yticks(fontname = "Times New Roman")  # This argument will change the font.
pl.xticks(fontname = "Times New Roman")  # This argument will change the font.
pl.rcParams['legend.fontsize'] = 18.0
pl.xlabel('triangle area / $nm^2$',fontsize=22, fontname = "Times New Roman")
pl.ylabel('frequency',fontsize=22, fontname = "Times New Roman")
pl.ylim(0,45)
pl.xlim(0,5)
pl.legend()
fig.tight_layout()

fig.savefig('/home/lukas/master_thesis_newdesign/figures/vdb_bins50_erboco_faces_1000dpi.eps', format='eps', dpi=1000)


pl.figure()
fig2, ax = pl.subplots()

pl.hist(ivas_triangle_areas, bins=bins,label='MC', color='blue')
pl.grid()
pl.tick_params(axis='both', which='major', labelsize=19)
pl.yticks(fontname = "Times New Roman")  # This argument will change the font.
pl.xticks(fontname = "Times New Roman")  # This argument will change the font.
pl.rcParams['legend.fontsize'] = 18.0
pl.xlim(0,5)
pl.ylim(0,45)
pl.xlabel('triangle area / $nm^2$',fontsize=22, fontname = "Times New Roman")
pl.ylabel('frequency',fontsize=22, fontname = "Times New Roman")
pl.legend()
fig2.tight_layout()

fig2.savefig('/home/lukas/master_thesis_newdesign/figures/ivas_bins50_erboco_faces_1000dpi.eps', format='eps', dpi=1000)

#pl.tick_params(axis='both', which='major', labelsize=19)
#pl.yticks(fontname = "Times New Roman")  # This argument will change the font.
#pl.xticks(fontname = "Times New Roman")  # This argument will change the font.
#pl.xlabel('triangle index ' ,fontsize=22, fontname = "Times New Roman")
#pl.ylabel('scalar triple product', fontsize=22, fontname = "Times New Roman")
#pl.grid()
#pl.rcParams['legend.fontsize'] = 18.0
#pl.legend()