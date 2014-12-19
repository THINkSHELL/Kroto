#!/usr/bin/env python
# -*- coding: utf-8 -*-

import clr
clr.AddReference("mtrand")

import numpy as np
import rhinoscriptsyntax as rs
import Rhino

def upwardFace(x1, x2, x3):
    v2 = rs.VectorCreate(x2, x1)
    v3 = rs.VectorCreate(x3, x1)
    vv = rs.VectorCrossProduct(v2, v3)
    z = rs.VectorCreate([0,0,1], [0,0,0])
    a = rs.VectorDotProduct(vv, z)
    return a > 0

def orientMeshFaces(mesh, connec):
    vertices = rs.MeshVertices(mesh)
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
    if debug: print vertices[i]
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
    newVertex = ( vertices[i]
                + speed * ((qMi.getI() * qMiX2).getT()-vertices[i]) )
    if (((newVertex-vertices[i])*(newVertex-vertices[i]).getT())[0, 0]
        > var):
        var = ((newVertex - vertices[i])*(newVertex - vertices[i])
               .getT())[0, 0]
    return var, newVertex
    
def iterateOneStep(vertices, vertexFaces, naked, q, iter, speed):
    newVertices = np.matrix(np.copy(vertices))
    iter += 1
    var = 0
    for i in range(len(vertices)):
        if not naked[i]:
            var, newVertices[i] = iterateVertex(i, vertices, vertexFaces,
                                                q, var, speed)
    vertices = np.matrix(np.copy(newVertices))
    return iter, var, vertices
    
def meshDistance(vertices, objective):
    distance = 0
    for i in range(len(vertices)):
        if (((objective[i]-vertices[i])*(objective[i]-vertices[i]).getT())[0, 0]
            > distance):
            distance = ((objective[i] - vertices[i])*(objective[i] - vertices[i])
                   .getT())[0, 0]
    return distance
    
def minimizeMesh(mesh, q=None, itermax=10, lmax=0.01, reference=None, speed=0.5):
    if not q: q = np.ones(rs.MeshFaceCount(mesh))
    if graphic or showResult: meshi = mesh
    oldVertices = rs.MeshVertices(mesh)
    vertices = np.matrix(np.zeros((3, len(oldVertices)))).getT()
    if not reference: reference = vertices
    
    for i, v in enumerate(oldVertices):
        vertices[i, 0] = v[0]
        vertices[i, 1] = v[1]
        vertices[i, 2] = v[2]
    connec = rs.MeshFaceVertices(mesh)
    vertexFaces = orientMeshFaces(mesh, connec)
    naked = rs.MeshNakedEdgePoints(mesh)
    
    iter = 0
    var = 2*lmax
    distances = []
    while (iter < itermax) & (var > lmax):
        iter, var, vertices = iterateOneStep(vertices, vertexFaces, naked,
                                             q, iter, speed)
        distances.append(meshDistance(vertices, reference))
        if graphic:
            rs.HideObject(meshi)
            V = []
            for i in range(0, len(vertices)):
                V.append((vertices[i, 0], vertices[i, 1], vertices[i, 2]))
            meshi = rs.AddMesh(V, connec)
            print (u"itération {},"
                   u"\nDéplacement^2 maximum depuis "
                   u"l'itération précedente : {} mm²").format(iter, var)
        else:
            print var
    
    if showResult:
        rs.HideObject(meshi)
        V = []
        for i in range(0, len(vertices)):
            V.append((vertices[i, 0], vertices[i, 1], vertices[i, 2]))
        rs.AddMesh(V, connec)
        
    return vertices, distances

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
