#!/usr/bin/env python
# -*- coding: utf-8 -*-

from imp import reload
import MeshMinimize as mm
import rhinoscriptsyntax as rs

reload(mm)

mm.debug = 0
mm.graphic = 0
mm.showResult = 1

lmax = 0.001

itermax = 5

mesh = rs.GetObject("select mesh", rs.filter.mesh)

objective, distances = mm.minimizeMesh(mesh, q=None, itermax=itermax,
                                    lmax=lmax, reference=None, speed=0.7)


    

"""
numspeeds = 10
distance = []

for i in range(numspeeds):
    speed = (i+1)/numspeeds
    convergenceStudy(speed)
"""

# convergenceStudy(0.75)
# convergenceStudy(0.76)

# minimizeMesh(vertices, vertexFaces, q, itermax, lmax, vertices, 0.5)
    
