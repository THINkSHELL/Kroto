#!/usr/bin/env python
# -*- coding: utf-8 -*-

from imp import reload
import MeshMinimize as mm
import rhinoscriptsyntax as rs
import timeit
import clr

reload(mm)

clr.EnableProfiler(True)

mm.debug = 0
mm.graphic = 0
mm.showResult = 0

lmax = 0.001
itermax = 1000

mesh = rs.GetObject("select mesh", rs.filter.mesh)

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped
    
grad = wrapper(mm.minimizeMesh, mesh, q=None, itermax=itermax,
                                    lmax=lmax, reference=None, speed=0.7,
                                    method='gradient')

fix = wrapper(mm.minimizeMesh, mesh, q=None, itermax=itermax,
                                    lmax=lmax, reference=None, speed=1,
                                    method='fixed-point')

def print_profile():
    for p in clr.GetProfilerData():
        print 'done'
        print '%s\t%d\t%d\t%d' % (p.Name, p.InclusiveTime, p.ExclusiveTime, p.Calls)

objective, distances = mm.minimizeMesh(mesh, q=None, itermax=itermax,
                                    lmax=lmax, reference=None, speed=1,
                                    method='fixed-point')

#print_profile()

"""
timegrad = timeit.timeit(grad, number=3)
timefix = timeit.timeit(fix, number=3)

print 'grad : {}'.format(timegrad)
print 'fix : {}'.format(timefix)
"""
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
    
