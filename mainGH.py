#!/usr/bin/env python
# -*- coding: utf-8 -*-

from imp import reload
import MeshMinimize as mm
import rhinoscriptsyntax as rs
import ghpythonlib.components as ghcomp

reload(mm)

print options

for key in options.__dict__:
    mm.__dict__[key] = options.__dict__[key]

qCables = None
if cables: qCables = [1000 for i in cables]

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped
    
call = wrapper(mm.minimizeMesh, mesh, cables, fixed=None, qs=None, qCables=qCables, reference=None)

if toggle:
    vertices = mm.minimizeMesh(mesh, cables, fixed=None, qs=None, qCables=qCables, reference=None)

meshOut = ghcomp.ConstructMesh(objective, rs.MeshFaceVertices(mesh))

displacements = [rs.VectorCreate(v1, v2) for (v1, v2) in zip(rs.MeshVertices(mesh), vertices)]