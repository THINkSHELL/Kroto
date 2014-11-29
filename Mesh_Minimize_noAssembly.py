#!/usr/bin/env python
# -*- coding: utf-8 -*-

import clr
clr.AddReference("mtrand")

import numpy as np
import rhinoscriptsyntax as rs
import Rhino

mesh = rs.GetObject("select mesh", rs.filter.mesh)
meshi = mesh

q = np.ones(rs.MeshFaceCount(mesh))

lmax = 1

connec = rs.MeshFaceVertices(mesh)
oldVertices = rs.MeshVertices(mesh)
naked = rs.MeshNakedEdgePoints(mesh)
vertexFaces = {}

vertices = np.matrix(np.zeros((3, len(oldVertices)))).getT()
for i, v in enumerate(oldVertices):
    vertices[i, 0] = v[0]
    vertices[i, 1] = v[1]
    vertices[i, 2] = v[2]
newVertices = np.matrix(np.copy(vertices))
var = 2*lmax
itermax = 50
iter = 0
speed = 1/2

def upwardFace(x1, x2, x3):
    v2 = rs.VectorCreate(x2, x1)
    v3 = rs.VectorCreate(x3, x1)
    vv = rs.VectorCrossProduct(v2, v3)
    z = rs.VectorCreate([0,0,1], [0,0,0])
    a = rs.VectorDotProduct(vv, z)
    return a > 0

for i, vertex in enumerate(oldVertices):
    for j, face in enumerate(rs.MeshVertexFaces(mesh, i)):
        others = list(set(connec[face]) - {i, })
        [x2, x3] = [oldVertices[n] for n in others]
        if upwardFace(vertex, x2, x3):
            vertexFaces[(i, j)] = [i, others[0], others[1]]
            # print 'not flip'
        else:
            vertexFaces[(i, j)] = [i, others[1], others[0]]
            # print 'flip'


while (iter < itermax) & (var > lmax):
    iter += 1
    var = 0
    for i, vertex in enumerate(vertices):
        if not naked[i]:
            qMi = np.matrix(np.zeros([3, 3]))
            qMiX2 = np.matrix(np.zeros([3, 1]))
            j = 0
            # print '####'
            # print vertex
            while vertexFaces.get((i, j), 0):
                x2 = vertices[vertexFaces[(i, j)][1]]
                x3 = vertices[vertexFaces[(i, j)][2]]
                x23 = x3 - x2
                qMij = q[face] * ( (x23 * x23.getT())[0, 0] 
                                    * np.matrix(np.eye(3))
                               - np.matrix([[x23[0, k]*x23[0, l] 
                                            for k in range (0,3)]
                                            for l in range(0, 3)]
                                          ))
                qMi += qMij
                qMiX2 += qMij * np.matrix([x2[0, 0], x2[0, 1], x2[0, 2]]).getT()
                # print (i, j)
                # print x2, x3
                # print qMij
                # print qMij * np.matrix([x2[0, 0], x2[0, 1], x2[0, 2]]).getT()
                # print '---'
                j += 1
            # print qMi
            # print qMiX2
            newVertices[i] = ( vertices[i]
                        + speed * ((qMi.getI() * qMiX2).getT()-vertices[i]) )
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
    print (u"itération {},"
           u"\nDéplacement^2 maximum depuis "
           u"l'itération précedente : {} mm²").format(iter, var)
    
# rs.HideObject(meshi)
# V = []
# for i in range(0, len(newVertices)):
#     V.append((newVertices[i, 0], newVertices[i, 1], newVertices[i, 2]))
# rs.AddMesh(V, connec)