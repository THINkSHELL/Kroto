#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import vectorWorks as vw
import rhinoscriptsyntax as rs
import Rhino

debug = 0
graphic = 0
showResult = 1
lmax = 0.01
itermax = 10
speed = 1
method = 'fixed-point'
fixedCableEnds = True


def upwardFace(x1, x2, x3):
    # ??? unnecessary, equations are already invariant by faces permutations
    # around a vertex ???
    """
    v2 = rs.VectorCreate(x2, x1)
    v3 = rs.VectorCreate(x3, x1)
    vv = rs.VectorCrossProduct(v2, v3)
    z = rs.VectorCreate([0,0,1], [0,0,0])
    a = rs.VectorDotProduct(vv, z)
    return a > 0
    """
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

def iterateVertex(i, vertices, vertexFaces, qs, ql, nCable, var):
    qMi = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    qMiX2 = [0, 0, 0]
    ql2i = 0
    qli = sum(ql[i])
    qliX2 = [0, 0, 0]
    j = 0
    if debug: print '####'
    if debug: print vertices[i]
    while vertexFaces.get((i, j), 0):
        x2 = vertices[vertexFaces[(i, j)][1]]
        x3 = vertices[vertexFaces[(i, j)][2]]
        x23 = vw.matminus(x3, x2)
        qij = qs[vertexFaces[(i, j)][1]]
        qMij = vw.matmul(qij, 
                  vw.matminus(
                    vw.matmul(vw.dotproduct(x23, x23), vw.id),
                    vw.veckronproduct(x23, x23)
                           )
                        )
        qMi = vw.matplus(qMi, qMij)
        qMiX2 = vw.matplus(qMiX2, vw.matmul(qMij, x2))
        ql2i += qij * vw.dotproduct(x23, x23)
        if debug: 
            print ''.join([str((i, j)), '\n', str(x2), str(x3), '\n', str(x23), '\n', str(ql), '\n', str(qMij)])
            print str(vw.matmul(qMij, x2)) + '\n' + str(ql2i) + '\n' + '---'
        j += 1
    for j, v in enumerate(nCable[i]):
        qliX2 = vw.matplus(qliX2, vw.matmul(ql[i][j], vertices[v]))
    if debug: print qMi + '\n' + qMiX2
    if method == 'gradient':
        newVertex = vw.matplus(
                        vertices[i],
                        vw.matmul(
                            speed,
                            vw.matminus(
                                vw.matmul(
                                    vw.inverse(
                                        vw.matplus(
                                            qMi,
                                            vw.matmul(qli, vw.id)
                                        )
                                    ),
                                    vw.matplus(
                                        qMiX2,
                                        qliX2
                                    )
                                ),
                                vertices[i]
                            )
                        )
                    )
    elif method == 'fixed-point':
        newVertex = vw.matplus( 
                        vw.matmul(
                            speed/(ql2i+qli), 
                            vw.matminus(
                                vw.matplus(
                                    qMiX2,
                                    qliX2
                                ),
                                vw.matmul(
                                    vw.matplus(
                                        qMi,
                                        vw.matmul(qli, vw.id)
                                    ),
                                    vertices[i]
                                )
                            )
                        ),
                        vertices[i]
                    )
    if (vw.dotproduct(vw.matminus(newVertex, vertices[i]), 
                       vw.matminus(newVertex, vertices[i]))
        > var):
        var = vw.dotproduct(vw.matminus(newVertex, vertices[i]), 
                       vw.matminus(newVertex, vertices[i]))
    return var, newVertex
    
def iterateOneStep(vertices, vertexFaces, fixed, qs, ql, nCable, iter):
    newVertices = copy.deepcopy(vertices)
    iter += 1
    var = 0
    for i in range(len(vertices)):
        if not fixed[i]:
            var, newVertices[i] = iterateVertex(i, vertices, vertexFaces, qs, ql, nCable, var)
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
    
def defineCables(cables, qCables, vertices, naked, fixed):
    tol = rs.UnitAbsoluteTolerance()
    temp = [[] for i in cables]
    ql = [[] for i in vertices]
    nCable = [[] for i in vertices]
    for v, vertex in enumerate(vertices):
        if naked[v]:
            for i, cable in enumerate(cables):
                if rs.Distance(
                   rs.EvaluateCurve(cable, rs.CurveClosestPoint(cable, vertex)),
                   vertex
                   ) < tol:
                    temp[i].append(v)
                    if ( fixedCableEnds and 
                         min(rs.Distance(rs.CurveEndPoint(cable), vertex),
                           rs.Distance(rs.CurveStartPoint(cable), vertex)) < tol
                       ):
                        if debug: rs.AddPoint(vertex)
                        fixed[v] = True
    for i, cable in enumerate(cables):
        temp[i].sort(key = lambda v: rs.CurveClosestPoint(cable, vertices[v]))
        for j in range(len(temp[i])):
            if ( (not fixed[temp[i][j]]) and j and (j-len(temp[i])+1) ):
                nCable[temp[i][j]].append(temp[i][j-1])
                nCable[temp[i][j]].append(temp[i][j+1])
                ql[temp[i][j]].append(qCables[i])
                ql[temp[i][j]].append(qCables[i])
    return ql, nCable, fixed
    
def minimizeMesh(mesh, cables=None, fixed=None, qs=None, qCables=None, reference=None):
    if not qs: qs = [1 for i in range(rs.MeshFaceCount(mesh))]
    if graphic or showResult: meshi = mesh
    oldVertices = rs.MeshVertices(mesh)
    vertices = [[oldVertices[i][0], oldVertices[i][1], oldVertices[i][2]] 
                for i in range(len(oldVertices))]
    if not reference: reference = vertices
    connec = rs.MeshFaceVertices(mesh)
    vertexFaces = orientMeshFaces(mesh, connec)
    naked = rs.MeshNakedEdgePoints(mesh)
    if not fixed:
        if not cables:
            fixed = naked
        else:
            fixed = [False for i in vertices]
    if cables:
        if not qCables: qCables = [1 for i in cables]
        ql, nCable, fixed = defineCables(cables, qCables, oldVertices, naked, fixed)
            
    iter = 0
    var = 2*lmax
    while (iter < itermax) & (var > lmax):
        iter, var, vertices = iterateOneStep(vertices, vertexFaces, fixed, qs, ql, nCable, iter)
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
        
    return vertices