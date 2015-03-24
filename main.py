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
mm.SAVE_RESULTS = 0
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
if cables:
    q_cables = [1000 for i in cables]

pressure = 100


def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped


def print_profile():
    for p in clr.GetProfilerData():
        print 'done'
        print '%s\t%d\t%d\t%d' % (p.Name, p.InclusiveTime, p.ExclusiveTime,
                                  p.Calls)

call = wrapper(mm.minimize_mesh, mesh, cables, fixed=None, qs=None,
               q_cables=q_cables, reference=None)

objective = mm.minimize_mesh(mesh, cables, fixed=None, qs=None,
                             q_cables=q_cables, reference=None, P=pressure)
