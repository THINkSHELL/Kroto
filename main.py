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
mm.showResult = 1
mm.lmax = 0.001
mm.itermax = 100
mm.speed = 1
mm.method = 'fixed-point'
mm.fixedCableEnds = True


mesh = rs.GetObject("Select mesh", rs.filter.mesh)
cables = rs.GetObjects("Select cables", rs.filter.curve)

qCables = None
if cables: qCables = [1000 for i in cables]

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped
    
call = wrapper(mm.minimizeMesh, mesh, cables, fixed=None, qs=None, qCables=qCables, reference=None)

def print_profile():
    for p in clr.GetProfilerData():
        print 'done'
        print '%s\t%d\t%d\t%d' % (p.Name, p.InclusiveTime, p.ExclusiveTime, p.Calls)

objective = mm.minimizeMesh(mesh, cables, fixed=None, qs=None, qCables=qCables, reference=None)

#print_profile()

"""
mm.speed = 0.7
mm.method = 'gradient'
timegrad = timeit.timeit(call, number=3)
mm.speed = 1
mm.method = 'fixed-point'
timefix = timeit.timeit(call, number=3)

print 'grad : {}'.format(timegrad)
print 'fix : {}'.format(timefix)
"""
"""

def convergenceStudy(speed):
    print '----'
    print speed
    distance = minimizeMesh(
                vertices, vertexFaces, q, 100, 0.001, vertices, speed)[1][::-1]
    print distance
    try:
        np.savetxt(
            "C:\\Users\\Pierre\\Desktop\\" + speed + ".csv",
            distance,
            delimiter = ",")
    except:
        pass


numspeeds = 10
distance = []

for i in range(numspeeds):
    speed = (i+1)/numspeeds
    convergenceStudy(speed)
"""

# convergenceStudy(0.75)
# convergenceStudy(0.76)

# minimizeMesh(vertices, vertexFaces, q, itermax, lmax, vertices, 0.5)
    
