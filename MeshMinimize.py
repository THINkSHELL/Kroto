#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Implements the Stress Density Method as a Rhino mesh solver
Takes a problem defined as a triangular mesh (working with quadrangular 
meshes as well would take a lot of rewrite), edges conditions and solver
options.

Returns a mesh close to a minimal surface. The result is a true minimal
surface if the edges are fixed, pressure null and the method iterated 
until density coefficients yield a uniform stress field.
"""

import copy
import vectorworks as vw
import rhinoscriptsyntax as rs
import meshminimizehelpers as mmh


# Define default options for the solver
# DEBUG = verbose option
# GRAPHIC = build every itereation mesh in Rhino (stand-alone script)
# SHOW_RESULT = build the resulting mesh in Rhino (stand-alone script)
# MAX_DISP = maximum displacement for the break criterion of the
#           convergence algorithm
# MAX_ITER = maximum number of iterations in the fixed-qs convergence 
#            algorithm
# MAX_ITER_QS = maximum number of iterations for the qs convergence 
#               algorithm
# MAX_DEV_SIGMA = maximum deviation around mean-value for the surface 
#                 stresses
# SPEED = displacement speed at each iteration (necessary for stability
#         of the gradient method)
# METHOD = 'fixed-point' or 'gradient'
# FIXED_CABLE_ENDS = cables ends are considered fixed
DEBUG = 0
GRAPHIC = 0
SHOW_RESULT = 0
MAX_DISP = 0.01
MAX_ITER = 10
MAX_ITER_QS = 10
MAX_DEV_SIGMA = 1
SPEED = 1
METHOD = 'fixed-point'
FIXED_CABLE_ENDS = True

def iterate_vertex(i, vertices, vertex_faces, naked, qs, ql, n_cable, P6, var):
    """Updates the position of a node in the mesh
    
    Arguments:
      i = the node number in the mesh
      vertices = the list of the mesh nodes
      vertex_faces = the list of faces adjacent to a node
      naked = list of booleans, True if vertex is on a naked edge
      qs = the stress density for all faces in the mesh
      ql = the cable force density
      n_cable = the list of list of cable segments from each vertex
      var = the current maximum displacement in the iteration
    Returns
      var = the updated maximum displacement
      new_vertex = the updated position of the node
    """
    
    """First we initialize the intermediate vectors for calculation
    All numberings (usually indexed by j) are relative to the current 
    vertex, i.
    
    .  = vector dot-product ( [n,1] . [n,1] -> [1,1] )
    /\ = vector cross-product ( [n,1] /\ [n,1] -> [n,1] )
    x  = vector kronecker product ( [n,1] x [n,1] -> [n,n] )
    *  = classical multipliation of scalars, vectors or matrices
    
    m_i = faces around the vertex i
    n_i = cables around the vertex i
    X_2j, X_3j = points 2 and 3 of the face #j around the vertex i, 
                 [3x1] vector
    X_2j3j = X_3j - X_2j = vector from 2j to 3j, [3x1] vector
    
    M_(i,j) = (X_2j3j . X_2j3j)*Id - (X_2j3j x X_2j3j)
              [3,3] matrix representing the dependant faces areas around
              the vertex i
    
    qs_j  = face number j surface stress density coefficient
    ql_j  = cable segment number j force density coefficient (for points
            in the middle of a cable, count each side once)
    P6     = pressure / 3 (uniform scalar at the moment)
    
    qMi   = sum(j = [1, m_i]; qs_j * M_(i,j))
          = local stiffness matrix, [3x3] matrix
    qMiX2 = sum(j = [1, m_i]; qs_j * M_(i,j) * X_2j) 
          = local membrane forces on the vertex, [3x1] vector
    ql2i  = sum(j = [1, m_i]; qs_j * l_ij^2) 
          = local membrane stifness coefficient, scalar
    qli   = sum(j = [1, n_i]; ql_j)
          = local cable stifness coefficient, scalar
    qliX2 = sum(j = [1, n_i]; ql_j * M_(i,j) * X_2j) 
          = local cable forces on the vertex, [3x1] vector
    PX2X3 = sum(j = [1, m_i]; P/6 * (X_2j/\X_3j + X_2j3j/\X_i)
          = pressure membrane forces around the vertex, [3x1] vector
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
    
    if DEBUG: print '####'
    if DEBUG: print vertices[i]
    
    # Iterate for faces j around current vertex i
    while vertex_faces.get((i, j), 0): # Face exists in vertex_faces table
        x2 = vertices[vertex_faces[(i, j)][1]] # = X_2j, 3D point, [3x1] vector
        x3 = vertices[vertex_faces[(i, j)][2]] # = X_3j, 3D point, [3x1] vector
        x23 = vw.matminus(x3, x2) # = vector X_2j3j, [3x1] vector
        qij = qs[vertex_faces[(i, j)][1]] # = qs_j = Surface stress density coef
        
        # qMij = qs_j * M_(i,j)
        # M_(i,j) = (X_2j3j . X_2j3j)*Id - (X_2j3j x X_2j3j)
        qMij = vw.matmul(qij, 
                  vw.matminus(
                    vw.matmul(vw.dotproduct(x23, x23), vw.id),
                    vw.veckronproduct(x23, x23)
                           )
                        )
                        
        # update qMi and qMiX2 running sums
        qMi = vw.matplus(qMi, qMij)
        qMiX2 = vw.matplus(qMiX2, vw.matmul(qMij, x2))
        
        # ql2i = sum(j = [1, m_i]; qs_j * l_ij^2), running sum update, scalar
        ql2i += qij * vw.dotproduct(x23, x23)
        
        # PX2X3 = sum(j = [1, m_i]; P/6 * (X_2j/\X_3j + X_2j3j/\X_i), 
        # running sum update, [3,1] vector
        PX2X3 = vw.matplus(PX2X3, vw.matmul(P6, vw.crossproduct(x2, x3)))
        # The contour term is non-zero only for edge nodes
        if naked[i]:
            PX2X3 = vw.matplus(PX2X3, 
                               vw.matmul(P6, vw.crossproduct(x23, vertices[i])))
        
        if DEBUG: 
            print '\n'.join(
                          [str((i, j)), str(x2), str(x3), str(x23), str(ql),
                          str(qMij), vw.matmul(qMij, x2), str(ql2i), '---'])
        j += 1
        #END while
        
    if DEBUG:
        print '\n'.join([qMi, qMiX2, P, PX2X3])
        
    # Iterate for cables segments j around vertex i
    if n_cable:
        for j, v in enumerate(n_cable[i]):
            # qliX2 = qliX2 = sum( j = [1, n_i]; ql_j * M_(i,j) * X_2j ),
            # Running sum update, [3x1] vector.
            qliX2 = vw.matplus(qliX2, vw.matmul(ql[i][j], vertices[v]))
            
    # Do the actual work here, both methods
    # X_i(t+1) = X_i(t) + SPEED / (qMi + qli) * 
    #                        ((qMiX2 + qliX2 + PX2X3) - X_i(t))
    if METHOD == 'gradient':
        new_vertex = vw.matplus(
                        vertices[i],
                        vw.matmul(
                            SPEED,
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
    # X_i(t+1) = X_i(t) + SPEED / (ql2i + qli) * 
    #                               (qMiX2 + qliX2 + PX2X3 - (qMi + qli).X_i(t))
    elif METHOD == 'fixed-point':
        new_vertex = vw.matplus(
                      vertices[i],
                        vw.matmul(
                            SPEED/(ql2i+qli), 
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
                    
    # Update displacement criterion, if needed.
    # max(squared disp) is not a really good criterion, to be improved..
    if (vw.dotproduct(vw.matminus(new_vertex, vertices[i]), 
                       vw.matminus(new_vertex, vertices[i]))
        > var):
        var = vw.dotproduct(vw.matminus(new_vertex, vertices[i]), 
                       vw.matminus(new_vertex, vertices[i]))
                       
    return var, new_vertex


def iterate_one_step(vertices, vertex_faces, naked, fixed, qs, ql, n_cable, P6,
                     iter):
    """Updates all nodes on the mesh once
    
    Arguments:
      vertices = list of mesh vertices, from rs.MeshVertices(mesh)
      vertex_faces = list of list of faces connected a vertex,
                     from orient_mesh_faces(mesh)
      naked = list of booleans, True if vertex is on a naked edge
      fixed = list of booleans, True if vertex is fixed
      qs = list of surface stress density coefficients for each face
      ql = list of cable force density, for each cable segment
      n_cable = list of list of cables segments connected to a vertex
      P6 = pressure / 3
      iter = current iteration number
    Returns:
      iter = updated iteration number
      var = maximum squared displacement for this iteration
      vertices = updated positions of vertices
    """
    
    # Copy vertices to a new list, so that we do not overwrite it
    new_vertices = copy.deepcopy(vertices)
    iter += 1
    var = 0
    for i in range(len(vertices)):
        if not fixed[i]:
            var, new_vertices[i] = iterate_vertex(i, vertices, vertex_faces, 
                                                 naked, qs, ql, n_cable, P6, 
                                                 var)
    vertices = copy.deepcopy(new_vertices)
    return iter, var, vertices


def iterate_fixed_qs(vertices, vertex_faces, connec, naked, fixed, qs, ql,
                     n_cable, P, meshi=None):
    """Iterates the problem with qs fixed, until a pseudo-minimal 
    surface is found.
    Arguments:
      vertices = list of mesh vertices, from rs.MeshVertices(mesh)
      vertex_faces = list of list of faces connected a vertex,
                     from orient_mesh_faces(mesh)
      connec = connectivity matrix in the mesh, from 
               rs.MeshFaceVertices(mesh)
      naked = list of booleans, True if vertex is on a naked edge
      fixed = list of booleans, True if vertex is fixed
      qs = list of surface stress density coefficients for each face
      ql = list of cable force density, for each cable segment
      n_cable = list of list of cable segments attached to each vertex
      P = pressure
      meshi = save-state for graphical display
    Returns:
      vertices = updated positions of vertices
    """
    
    # Initialize loop
    iter = 0
    var = 2*MAX_DISP
    
    # Loop while we can, get some display if wanted
    while (iter < MAX_ITER) & (var > MAX_DISP):
        iter, var, vertices = iterate_one_step(vertices, vertex_faces, 
                                               naked, fixed, qs, ql,
                                               n_cable, P/6, iter)
        if GRAPHIC:
            rs.HideObject(meshi)
            meshi = rs.AddMesh(vertices, connec)
            print (u"itération {},"
                   u"\nDéplacement^2 maximum depuis "
                   u"l'itération précedente : {} mm²").format(iter, var)
        else:
            print var
    
    return vertices


def minimize_mesh(mesh, cables=None, fixed=None, qs=None, q_cables=None,
                  reference=None, ql=None, n_cable=None, P=0):
    """Iterates a mesh until it is close to a minimal surface.
    
    Arguments:
      mesh = mesh to calculate, Rhino GUID
      cables = polylines representing the cables
      fixed = list of booleans, True if vertex is fixed
      qs = list of surface stress density coefficients for each face
      q_cables = list of force density coefficients for each cable
      reference = reference mesh for comparisons, unused
      ql = list of cable force density, for each cable segment
      n_cable = list of list of cable segments connected to a vertex
      P = pressure
    Returns:
        vertices = vertices at new position
    """
    
    #
    # Initialize
    # qs defaults to 1 everywhere
    # meshi remembers current mesh if we need to hide it in Rhino
    # vertices strips RhinoCommon's 3D points to a bare 3-vector
    # reference defaults to the initial mesh
    # vertex_faces is reordered from Rhino's mesh
    # naked is extracted from Rhino
    # fixed defaults to naked if edges are fixed, otherwise defined in
    #                   mmh.define_cables
    # ql, n_cable are defined by mmh.define_cables if necessary
    # q_cables defaults to 1
    #
    
    if not qs:
        qs = [1 for i in range(rs.MeshFaceCount(mesh))]
    if GRAPHIC or SHOW_RESULT:
        meshi = mesh
    old_vertices = rs.MeshVertices(mesh)
    vertices = [[old_vertices[i][0], old_vertices[i][1], old_vertices[i][2]] 
                for i in range(len(old_vertices))]
    if not reference:
        reference = vertices
    vertex_faces = mmh.orient_mesh_faces(mesh)
    naked = rs.MeshNakedEdgePoints(mesh)
    connec = rs.MeshFaceVertices(mesh)
    if not fixed:
        if not cables:
            fixed = naked
        else:
            fixed = [False for i in vertices]
    ql = None
    n_cable = None
    if cables:
        if not q_cables:
            q_cables = [1 for i in cables]
        if not (ql and n_cable and fixed):
            ql, n_cable, fixed = mmh.define_cables(cables, q_cables, 
                                                   old_vertices, naked, fixed)
    
    iter_qs = 0
    dev_sigma = 2*MAX_DEV_SIGMA
    
    while (iter_qs < MAX_ITER_QS) & (dev_sigma > MAX_DEV_SIGMA):
        vertices = iterate_fixed_qs(vertices, vertex_faces, connec, naked, 
                                    fixed, qs, ql, n_cable, P)
        iter_qs += 1
        dev_sigma, qs = mmh.update_qs(mesh, qs)
    
    if SHOW_RESULT:
        rs.HideObject(meshi)
        rs.AddMesh(vertices, connec)
        
    return vertices