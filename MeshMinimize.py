#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implements the Stress Density Method as a Rhino mesh solver
Takes a problem defined as a mesh, edges conditions and solver options
Returns a mesh close to a minimal surface
"""

import copy
import vectorWorks as vw
import rhinoscriptsyntax as rs

#Define default options for the solver
debug = 0
graphic = 0
showResult = 1
lmax = 0.01
itermax = 10
speed = 1
method = 'fixed-point'
fixedCableEnds = True

def upwardFace(n, x1, x2, x3):
    """
    Returns True if the face orientation defined by the order of the passed in 
    points is the same as that of the mesh, passed by the normal n
    """
    v2 = rs.VectorCreate(x2, x1)
    v3 = rs.VectorCreate(x3, x1)
    vv = rs.VectorCrossProduct(v2, v3)
    a = rs.VectorDotProduct(vv, n)
    return a > 0

def orientMeshFaces(mesh, connec):
    """
    Orients faces around a node in a mesh to a consistant order and normal direction
    Arguments:
        mesh the mesh in RhinoCommon type
        connec lists of connected points for each point in the mesh
    Returns:
        vertexFaces list of faces adjacent to a node
    """
    
    normals = rs.MeshFaceNormals(mesh)
    vertices = rs.MeshVertices(mesh)
    vertexFaces = {}
    for i, vertex in enumerate(vertices):
        for j, face in enumerate(rs.MeshVertexFaces(mesh, i)):
            #Find the three distinct points, 4th is redundant in triangular meshes
            #and remove current point from the set
            others = list(set(connec[face]) - {i, })
            [x2, x3] = [vertices[n] for n in others]
            if upwardFace(normals[i], vertex, x2, x3):
                vertexFaces[(i, j)] = [i, others[0], others[1]]
                if debug: print 'not flip'
            else:
                vertexFaces[(i, j)] = [i, others[1], others[0]]
                if debug: print 'flip'
    return vertexFaces

def iterateVertex(i, vertices, vertexFaces, qs, ql, nCable, P, var):
    """
    Update the position of a node in the mesh
    Arguments:
        i the node number in the mesh
        vertices the list of the mesh nodes
        vertexFaces the list of faces adjacent to a node
        qs the stress density for all faces in the mesh
        ql the cable force density
        nCable the list of cables
        var the current maximum displacement in the iteration
    Returns
        var: the updated maximum displacement
        newVertex: the updated position of the node
    """
    
    qMi = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    qMiX2 = [0, 0, 0]
    ql2i = 0
    if ql:
        qli = sum(ql[i])
    else:
        qli = 0
    qliX2 = [0, 0, 0]
    PX2X3 = [0, 0, 0]
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
        PX2X3 = vw.matplus(PX2X3, vw.matmul(P, vw.crossproduct(x2, x3)))
        if debug: 
            print ''.join([str((i, j)), '\n', str(x2), str(x3), '\n', str(x23), '\n', str(ql), '\n', str(qMij)])
            print str(vw.matmul(qMij, x2)) + '\n' + str(ql2i) + '\n' + '---'
        j += 1
    if nCable:
        for j, v in enumerate(nCable[i]):
            qliX2 = vw.matplus(qliX2, vw.matmul(ql[i][j], vertices[v]))
    if debug: print qMi + '\n' + qMiX2
    print P
    print PX2X3
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
                                        vw.matplus(qliX2, PX2X3)
                                    )
                                ),
                                vertices[i]
                            )
                        )
                    )
    elif method == 'fixed-point':
        newVertex = vw.matplus(
                      vertices[i],
                        vw.matmul(
                            speed/(ql2i+qli), 
                            vw.matminus(
                                vw.matplus(
                                    qMiX2,
                                    vw.matplus(qliX2, PX2X3)
                                ),
                                vw.matmul(
                                    vw.matplus(
                                        qMi,
                                        vw.matmul(qli, vw.id)
                                    ),
                                    vertices[i]
                                )
                            )
                        )
                    )
    if (vw.dotproduct(vw.matminus(newVertex, vertices[i]), 
                       vw.matminus(newVertex, vertices[i]))
        > var):
        var = vw.dotproduct(vw.matminus(newVertex, vertices[i]), 
                       vw.matminus(newVertex, vertices[i]))
    return var, newVertex
    
def iterateOneStep(vertices, vertexFaces, fixed, qs, ql, nCable, P, iter):
    newVertices = copy.deepcopy(vertices)
    iter += 1
    var = 0
    for i in range(len(vertices)):
        if not fixed[i]:
            var, newVertices[i] = iterateVertex(i, vertices, vertexFaces, qs, ql, nCable, P, var)
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
    
def minimizeMesh(mesh, cables=None, fixed=None, qs=None, qCables=None, reference=None, ql=None, nCable=None, P=0):
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
    ql = None
    nCable = None
    if cables:
        if not qCables: qCables = [1 for i in cables]
        if not (ql and nCable and fixed):
            ql, nCable, fixed = defineCables(cables, qCables, oldVertices, naked, fixed)
    
    iter = 0
    var = 2*lmax
    while (iter < itermax) & (var > lmax):
        iter, var, vertices = iterateOneStep(vertices, vertexFaces, fixed, qs, ql, nCable, P, iter)
        if graphic:
            rs.HideObject(meshi)
            meshi = rs.AddMesh(vertices, connec)
            print (u"itération {},"
                   u"\nDéplacement^2 maximum depuis "
                   u"l'itération précedente : {} mm²").format(iter, var)
        else:
            print var
    
    #if showResult:
    #    rs.HideObject(meshi)
    #    rs.AddMesh(vertices, connec)
        
    return vertices