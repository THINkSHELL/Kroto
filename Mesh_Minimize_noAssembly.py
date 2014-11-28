#!/usr/bin/env python
# -*- coding: utf-8 -*-

import clr
clr.AddReference("mtrand")

import numpy as np
import rhinoscriptsyntax as rs

mesh = rs.GetObject("select mesh", rs.filter.mesh)
meshi = mesh

q = np.ones(rs.MeshFaceCount(mesh))

lmax = 1

connec = rs.MeshFaceVertices(mesh)
oldVertices = rs.MeshVertices(mesh)
naked = rs.MeshNakedEdgePoints(mesh)

vertices = np.matrix(np.zeros((3, len(oldVertices)))).getT()
for i, v in enumerate(oldVertices):
    vertices[i, 0] = v[0]
    vertices[i, 1] = v[1]
    vertices[i, 2] = v[2]
newVertices = np.matrix(np.copy(vertices))
var = 2*lmax
itermax = 10
iter = 0

while (iter < itermax) & (var > lmax):
    iter += 1
    print (u"itération {},"
            + u"\nDéplacement^^2 maximum depuis "
            + u"l'itération précedente : {} mm²").format(iter, var)
    var = 0
    for i, vertex in enumerate(vertices):
        if not naked[i]:
            qMi = np.matrix(np.zeros([3, 3]))
            qMiX2 = np.matrix(np.zeros([3, 1]))

            for face in rs.MeshVertexFaces(mesh, i):
                others = list(set(connec[face]) - {i, })
                [x2, x3] = [vertices[n] for n in others]
                x23 = x3 - x2
                #PB de l'ordre x2/x3 ? Orientation du maillage par Rhino...
                qMij = q[face] * ( (x23 * x23.getT())[0, 0] 
                                    * np.matrix(np.eye(3))
                               - np.matrix([[x23[0, k]*x23[0, l] 
                                            for k in range (0,3)]
                                            for l in range(0, 3)]
                                          ))
                qMi += qMij
                qMiX2 += qMij * np.matrix([x2[0, 0], x2[0, 1], x2[0, 2]]).getT()
            newVertices[i] = (qMi.getI() * qMiX2).getT()
            if (((newVertices[i] - vertex)*(newVertices[i]-vertex).getT())[0, 0]
                > var):
                var = ((newVertices[i] - vertex)*(newVertices[i] - vertex)
                        .getT())[0, 0]
    vertices = np.matrix(np.copy(newVertices))
    rs.HideObject(meshi)
    V = []
    for i in range(0, len(newVertices)):
        V.append((newVertices[i, 0], newVertices[i, 1], newVertices[i, 2]))
    meshi = rs.AddMesh(V, connec)
    
rs.HideObject(meshi)
V = []
for i in range(0, len(newVertices)):
    V.append((newVertices[i, 0], newVertices[i, 1], newVertices[i, 2]))
rs.AddMesh(V, connec)