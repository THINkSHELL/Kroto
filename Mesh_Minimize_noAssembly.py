#!/usr/bin/env python
# -*- coding: utf-8 -*-

import clr
clr.AddReference("mtrand")

import numpy as np
import rhinoscriptsyntax as rs
import Rhino

debug = 0
graphic = 1
showResult = 0

lmax = 0.001

itermax = 5

mesh = rs.GetObject("select mesh", rs.filter.mesh)

q = np.ones(rs.MeshFaceCount(mesh))

if graphic: meshi = mesh

connec = rs.MeshFaceVertices(mesh)
oldVertices = rs.MeshVertices(mesh)
naked = rs.MeshNakedEdgePoints(mesh)

vertices = np.matrix(np.zeros((3, len(oldVertices)))).getT()
for i, v in enumerate(oldVertices):
    vertices[i, 0] = v[0]
    vertices[i, 1] = v[1]
    vertices[i, 2] = v[2]
newVertices = np.matrix(np.copy(vertices))

def upwardFace(x1, x2, x3):
    v2 = rs.VectorCreate(x2, x1)
    v3 = rs.VectorCreate(x3, x1)
    vv = rs.VectorCrossProduct(v2, v3)
    z = rs.VectorCreate([0,0,1], [0,0,0])
    a = rs.VectorDotProduct(vv, z)
    return a > 0

def orientMeshFaces(mesh):
    vertices = rs.MeshVertices(mesh)
    connec = rs.MeshFaceVertices(mesh)
    vertexFaces = {}
    for i, vertex in enumerate(vertices):
        for j, face in enumerate(rs.MeshVertexFaces(mesh, i)):
            others = list(set(connec[face]) - {i, })
            [x2, x3] = [vertices[n] for n in others]
            if upwardFace(vertex, x2, x3):
                vertexFaces[(i, j)] = [i, others[0], others[1]]
                if debug: print 'not flip'
            else:
                vertexFaces[(i, j)] = [i, others[1], others[0]]
                if debug: print 'flip'
    return vertexFaces

def iterateVertex(i, vertices, vertexFaces, q, var, speed):
    qMi = np.matrix(np.zeros([3, 3]))
    qMiX2 = np.matrix(np.zeros([3, 1]))
    j = 0
    if debug: print '####'
    if debug: print vertex
    while vertexFaces.get((i, j), 0):
        x2 = vertices[vertexFaces[(i, j)][1]]
        x3 = vertices[vertexFaces[(i, j)][2]]
        x23 = x3 - x2
        qMij = q[j] * ( (x23 * x23.getT())[0, 0] * np.matrix(np.eye(3))
                         - np.matrix([[x23[0, k]*x23[0, l] for k in range (0,3)]
                                                           for l in range(0, 3)]
                                  )
                          )
        qMi += qMij
        qMiX2 += qMij * np.matrix([x2[0, 0], x2[0, 1], x2[0, 2]]).getT()
        if debug: print (i, j)
        if debug: print x2, x3
        if debug: print qMij
        if debug: print qMij * np.matrix([x2[0, 0], x2[0, 1], x2[0, 2]]).getT()
        if debug: print '---'
        j += 1
    if debug: print qMi
    if debug: print qMiX2
    newVertices[i] = ( vertices[i]
                + speed * ((qMi.getI() * qMiX2).getT()-vertices[i]) )
    if (((newVertices[i]-vertices[i])*(newVertices[i]-vertices[i]).getT())[0, 0]
        > var):
        var = ((newVertices[i] - vertices[i])*(newVertices[i] - vertices[i])
               .getT())[0, 0]
    return var
    
def iterateOneStep(vertices, vertexFaces, q, iter, speed):
    global meshi
    iter += 1
    var = 0
    for i in range(len(vertices)):
        if not naked[i]:
            var = iterateVertex(i, vertices, vertexFaces, q, var, speed)
    vertices = np.matrix(np.copy(newVertices))
    if graphic:
        rs.HideObject(meshi)
        V = []
        for i in range(0, len(newVertices)):
            V.append((newVertices[i, 0], newVertices[i, 1], newVertices[i, 2]))
        meshi = rs.AddMesh(V, connec)
        print (u"itération {},"
               u"\nDéplacement^2 maximum depuis "
               u"l'itération précedente : {} mm²").format(iter, var)
    else:
        print var
    return iter, var, vertices
    
def meshDistance(vertices, objective):
    distance = 0
    for i in range(len(vertices)):
        if (((objective[i]-vertices[i])*(objective[i]-vertices[i]).getT())[0, 0]
            > distance):
            distance = ((objective[i] - vertices[i])*(objective[i] - vertices[i])
                   .getT())[0, 0]
    return distance
    
def minimizeMesh(vertices, vertexFaces, q, itermax, lmax, reference, speed):
    iter = 0
    var = 2*lmax
    distances = []
    while (iter < itermax) & (var > lmax):
        iter, var, vertices = iterateOneStep(vertices, vertexFaces, q, iter, speed)
        distances.append(meshDistance(vertices, reference))
    return vertices, distances
    
vertexFaces = orientMeshFaces(mesh)

objective, distances = minimizeMesh(vertices, vertexFaces, q, 100, 0.001, vertices, 0.5)
numspeeds = 10
distance = []

def convergenceStudy(speed):
    print '----'
    print speed
    distance = minimizeMesh(
                vertices, vertexFaces, q, 100, 0.001, vertices, speed)[1][::-1]
    print distance
    try:
        np.savetxt(
            "C:\\Users\\Pierre\\Desktop\\" + speed + ".csv",
            np.asarray(distance),
            delimiter = ",")
    except:
        pass

"""
for i in range(numspeeds):
    speed = (i+1)/numspeeds
    convergenceStudy(speed)
"""

# convergenceStudy(0.75)
# convergenceStudy(0.76)

# minimizeMesh(vertices, vertexFaces, q, itermax, lmax, vertices, 0.5)
    
if showResult:
    rs.HideObject(meshi)
    V = []
    for i in range(0, len(newVertices)):
        V.append((newVertices[i, 0], newVertices[i, 1], newVertices[i, 2]))
    rs.AddMesh(V, connec)