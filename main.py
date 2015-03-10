#!/usr/bin/env python
# -*- coding: utf-8 -*-

from imp import reload
import meshminimize as mm
import rhinoscriptsyntax as rs
import clr

reload(mm)

clr.EnableProfiler(True)

mm.DEBUG = 0
mm.GRAPHIC = 0
mm.SHOW_RESULT = 1
mm.MAX_DISP = 0.001
mm.MAX_ITER = 100
mm.MAX_ITER_QS = 10
mm.MAX_DEV_SIGMA = 1
mm.SPEED = 1
mm.METHOD = 'fixed-point'
mm.FIXED_CABLE_ENDS = True


mesh = rs.GetObject("Select mesh", rs.filter.mesh)
cables = rs.GetObjects("Select cables", rs.filter.curve)

q_cables = None
if cables: q_cables = [1000 for i in cables]

pressure = 100

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped
    
call = wrapper(mm.minimize_mesh, mesh, cables, fixed=None, qs=None, q_cables=q_cables, reference=None)

def print_profile():
    for p in clr.GetProfilerData():
        print 'done'
        print '%s\t%d\t%d\t%d' % (p.Name, p.InclusiveTime, p.ExclusiveTime, p.Calls)

objective = mm.minimize_mesh(mesh, cables, fixed=None, qs=None, q_cables=q_cables, reference=None, P=pressure)

#print_profile()

"""
mm.SPEED = 0.7
mm.METHOD = 'gradient'
timegrad = timeit.timeit(call, number=3)
mm.SPEED = 1
mm.METHOD = 'fixed-point'
timefix = timeit.timeit(call, number=3)

print 'grad : {}'.format(timegrad)
print 'fix : {}'.format(timefix)
"""
"""

def convergenceStudy(SPEED):
    print '----'
    print SPEED
    distance = minimize_mesh(
                vertices, vertexFaces, q, 100, 0.001, vertices, SPEED)[1][::-1]
    print distance
    try:
        np.savetxt(
            "C:\\Users\\Pierre\\Desktop\\" + SPEED + ".csv",
            distance,
            delimiter = ",")
    except:
        pass


numspeeds = 10
distance = []

for i in range(numspeeds):
    SPEED = (i+1)/numspeeds
    convergenceStudy(SPEED)
"""

# convergenceStudy(0.75)
# convergenceStudy(0.76)

# minimize_mesh(vertices, vertexFaces, q, MAX_ITER, MAX_DISP, vertices, 0.5)
    
