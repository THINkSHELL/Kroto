#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import vectorWorks as vw
import rhinoscriptsyntax as rs
import Rhino

def upwardFace(x1, x2, x3):
    """
    v2 = rs.VectorCreate(x2, x1)
    v3 = rs.VectorCreate(x3, x1)
    vv = rs.VectorCrossProduct(v2, v3)
    z = rs.VectorCreate([0,0,1], [0,0,0])
    a = rs.VectorDotProduct(vv, z)
    return a > 0
    """
    # unnecessary, eqautions are already invariant by faces permutations around 
    # a vertex
    return True

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

def iterateVertex(i, vertices, vertexFaces, q, var, speed, method):
    qMi = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    qMiX2 = [0, 0, 0]
    qli = 0
    j = 0
    if debug: print '####'
    if debug: print vertices[i]
    while vertexFaces.get((i, j), 0):
        x2 = vertices[vertexFaces[(i, j)][1]]
        x3 = vertices[vertexFaces[(i, j)][2]]
        x23 = vw.matminus(x3, x2)
        qij = q[vertexFaces[(i, j)][1]]
        qMij = vw.matmul(qij, 
                  vw.matminus(
                    vw.matmul(vw.dotproduct(x23, x23), vw.id),
                    vw.veckronproduct(x23, x23)
                           )
                        )
        qMi = vw.matplus(qMi, qMij)
        qMiX2 = vw.matplus(qMiX2, vw.matmul(qMij, x2))
        qli += qij * vw.dotproduct(x23, x23)
        if debug: 
            print (i, j)
            print x2, x3
            print x23
            print q
            print qMij
            print vw.matmul(qMij, x2)
            print qli
            print '---'
        j += 1
    if debug: print qMi
    if debug: print qMiX2
    if method == 'gradient':
        newVertex = vw.matplus(vertices[i], 
                        vw.matmul(speed, 
                            vw.matminus(
                                        vw.matmul(vw.inverse(qMi), qMiX2),
                                        vertices[i]
                                       )
                            ))
    elif method == 'fixed-point':
        newVertex = vw.matplus( vw.matmul(speed/qli, vw.matminus(qMiX2,
                                                vw.matmul(qMi, vertices[i]))),
                                vertices[i] )
    if (vw.dotproduct(vw.matminus(newVertex, vertices[i]), 
                       vw.matminus(newVertex, vertices[i]))
        > var):
        var = vw.dotproduct(vw.matminus(newVertex, vertices[i]), 
                       vw.matminus(newVertex, vertices[i]))
    return var, newVertex
    
def iterateOneStep(vertices, vertexFaces, naked, q, iter, speed, method):
    newVertices = copy.deepcopy(vertices)
    iter += 1
    var = 0
    for i in range(len(vertices)):
        if not naked[i]:
            var, newVertices[i] = iterateVertex(i, vertices, vertexFaces,
                                                q, var, speed, method)
    vertices = copy.deepcopy(newVertices)
    return iter, var, vertices
    
def meshDistance(vertices, objective):
    distance = 0
    for i in range(len(vertices)):
        if (vw.dotproduct(vw.matminus(objective[i], vertices[i]),
                          vw.matminus(objective[i], vertices[i]))
            > distance):
            distance = vw.dotproduct(vw.matminus(objective[i], vertices[i]),
                                     vw.matminus(objective[i], vertices[i]))
    return distance
    
def minimizeMesh(mesh, q=None, itermax=10, lmax=0.01, reference=None, speed=0.5,
                 method='gradient'):
    if not q: q = [1 for i in range(rs.MeshFaceCount(mesh))]
    if graphic or showResult: meshi = mesh
    oldVertices = rs.MeshVertices(mesh)
    vertices = [[oldVertices[i][0], oldVertices[i][1], oldVertices[i][2]] 
                for i in range(len(oldVertices))]
    if not reference: reference = vertices
    connec = rs.MeshFaceVertices(mesh)
    vertexFaces = orientMeshFaces(mesh, connec)
    naked = rs.MeshNakedEdgePoints(mesh)
    
    iter = 0
    var = 2*lmax
    distances = []
    while (iter < itermax) & (var > lmax):
        iter, var, vertices = iterateOneStep(vertices, vertexFaces, naked,
                                             q, iter, speed, method)
        #distances.append(meshDistance(vertices, reference))
        if graphic:
            rs.HideObject(meshi)
            meshi = rs.AddMesh(vertices, connec)
            print (u"itération {},"
                   u"\nDéplacement^2 maximum depuis "
                   u"l'itération précedente : {} mm²").format(iter, var)
        else:
            print var
    
    if showResult:
        rs.HideObject(meshi)
        rs.AddMesh(vertices, connec)
        
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
