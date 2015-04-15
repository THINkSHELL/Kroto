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

from vectorworks import *  # noqa
from Grasshopper.Kernel.Data import GH_Path
from Grasshopper import DataTree
import copy
import rhinoscriptsyntax as rs
import meshminimizehelpers as mmh


# Define default options for the solver
# DEBUG = verbose option
# GRAPHIC = build every iteration mesh in Rhino (stand-alone script)
# SHOW_RESULT = build the resulting mesh in Rhino (stand-alone script)
# SAVE_RESULTS = save the values of all the forces at each iteration
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
# METHOD = 'fixed-point', 'seidel' or 'gradient'
# FIXED_CABLE_ENDS = cables ends are considered fixed
DEBUG = 0
GRAPHIC = 0
SHOW_RESULT = 0
SAVE_RESULTS = 0
MAX_DISP = 0.01
MAX_ITER = 10
MAX_ITER_QS = 10
MAX_DEV_SIGMA = 1
SPEED = 1
METHOD = 'fixed-point'
FIXED_CABLE_ENDS = True


def iterate_vertex(i, vertices, vertices_faces_nodes, vertex_faces,  # noqa
                   naked, qs, ql, n_cable, P6, G6, var, iter_qs, iter, out):
    """Updates the position of a node in the mesh

    Arguments:
      i = the node number in the mesh
      vertices = the list of the mesh nodes
      vertices_faces_nodes = list of lists of list of nodes consituting
            the faces connected to a vertex, from orient_mesh_faces(mesh)
      vertex_faces = the list of faces adjacent to the node
      naked = True if the vertex is on a naked edge
      qs = the stress density for all faces in the mesh
      ql = the cable force density
      n_cable = the list of list of cable segments from each vertex
      P6 = pressure / 6
      G6 = surfacic weight density / 6
      var = the current maximum displacement in the iteration
      iter_qs = current iteration in the qs loop
      iter = current iteration in the stresses loop
      out = dictionnary of lists saving the forces acting on each
            vertex (qMiX2, qliX2 and PX2X3)
    Returns
      var = the updated maximum displacement
      new_vertex = the updated position of the node
      out = results saved with new vertex
    """

    # All numberings (usually indexed by j) are relative to the current
    # vertex, i.

    # .  = vector dot-product ( [n,1] . [n,1] -> [1,1] )
    # /\ = vector cross-product ( [n,1] /\ [n,1] -> [n,1] )
    # x  = vector kronecker product ( [n,1] x [n,1] -> [n,n] )
    # *  = classical multipliation of scalars, vectors or matrices

    # m_i = faces around the vertex i
    # n_i = cables around the vertex i
    # X_2j, X_3j = points 2 and 3 of the face #j around the vertex i,
    #              [3x1] vector
    # X_2j3j = X_3j - X_2j = vector from 2j to 3j, [3x1] vector

    # M_(i,j) = (X_2j3j . X_2j3j)*Id - (X_2j3j x X_2j3j)
    #           [3,3] matrix representing the dependant faces areas around
    #           the vertex i

    # qs_j  = face number j surface stress density coefficient
    # ql_j  = cable segment number j force density coefficient (for points
    #         in the middle of a cable, count each side once)
    # P6    = pressure / 6 (uniform scalar at the moment)
    # G6    = surfacic weight density / 6 (uniform scalar at the moment)

    # qMi   = sum(j = [1, m_i]; qs_j * M_(i,j))
    #       = local stiffness matrix, [3x3] matrix
    # qMiX2 = sum(j = [1, m_i]; qs_j * M_(i,j) * X_2j)
    #       = local membrane forces on the vertex, [3x1] vector
    # ql2i  = sum(j = [1, m_i]; qs_j * l_ij^2)
    #       = local membrane stifness coefficient, scalar
    # qli   = sum(j = [1, n_i]; ql_j)
    #       = local cable stifness coefficient, scalar
    # qliX2 = sum(j = [1, n_i]; ql_j * M_(i,j) * X_2j)
    #       = local cable forces on the vertex, [3x1] vector
    # PX2X3 = sum(j = [1, m_i]; P/6 * (X_2j/\X_3j + X_2j3j/\X_i))
    #       = pressure membrane forces around the vertex, [3x1] vector
    # G3    = sum(j = [1, m_i]; norm(X_12j/\X_13j))
    #       = u3 term / G6 of the gravity forces around the vertex, scalar
    # Fg    = [0, 0, - G3 * G6] = gravity forces, [3x1] vector

    # First we initialize the intermediate vectors for calculation
    qMi = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]  # noqa
    qMiX2 = [0, 0, 0]  # noqa
    ql2i = 0
    if ql:
        qli = sum(ql[i])
    else:
        qli = 0
    qliX2 = [0, 0, 0]  # noqa
    PX2X3 = [0, 0, 0]  # noqa
    G3 = 0  # noqa

    j = 0

    if DEBUG:
        print '####'
    if DEBUG:
        print vertices[i]

    # Iterate for faces j around current vertex i
    while vertices_faces_nodes.get((i, j), 0):  # Face is in vertex_faces table
        x2 = vertices[vertices_faces_nodes[(i, j)][1]]  # = X_2j,3D point,[3x1]
        x3 = vertices[vertices_faces_nodes[(i, j)][2]]  # = X_3j,3D point,[3x1]
        x23 = matminus(x3, x2)  # = vector X_2j3j, [3x1] vector
        if G6:
            x12 = matminus(x2, vertices[i])
            x13 = matminus(x3, vertices[i])
        qij = qs[vertex_faces[j]]  # = qs_j = Surface stress density coef

        # qMij = qs_j * M_(i,j)
        # M_(i,j) = (X_2j3j . X_2j3j)*Id - (X_2j3j x X_2j3j)
        qMij = matmul(  # noqa
            qij,
            matminus(
                diag(dotproduct(x23, x23), dim=3),
                veckronproduct(x23, x23)
            )
        )

        # update qMi and qMiX2 running sums
        qMi = matplus(qMi, qMij)  # noqa
        qMiX2 = matplus(qMiX2, matmul(qMij, x2))  # noqa

        # ql2i = sum(j = [1, m_i]; qs_j * l_ij^2), running sum update, scalar
        ql2i += qij * dotproduct(x23, x23)

        # PX2X3 = sum(j = [1, m_i]; P/6 * (X_2j/\X_3j + X_2j3j/\X_i),
        # running sum update, [3,1] vector
        if P6:
            PX2X3 = matplus(  # noqa
                PX2X3,  # noqa
                matmul(P6, crossproduct(x2, x3))
            )
            # The contour term is non-zero only for edge nodes
            if naked:
                PX2X3 = matplus(  # noqa
                    PX2X3,
                    matmul(P6, crossproduct(x23, vertices[i]))
                )

        # G3 = sum(j = [1, m_i]; norm(X_12j/\X_13j)),
        # running sum update, scalar
        if G6:
            G3 += dotproduct(
                crossproduct(x12, x13),
                crossproduct(x12, x13)
            ) ** 0.5

        if DEBUG:
            print '\n'.join([
                str((i, j)), str(x2), str(x3), str(x23), str(ql),
                str(qMij), matmul(qMij, x2), str(ql2i), '---'
            ])
        j += 1
        # END while

    if DEBUG:
        print '\n'.join([qMi, qMiX2, PX2X3])

    # Iterate for cables segments j around vertex i
    if n_cable:
        for j, v in enumerate(n_cable[i]):
            # qliX2 = qliX2 = sum( j = [1, n_i]; ql_j * M_(i,j) * X_2j ),
            # Running sum update, [3x1] vector.
            qliX2 = matplus(qliX2, matmul(ql[i][j], vertices[v]))  # noqa

    # Build the gravity forces vector
    # Fg = [0, 0, - G3 * G6]
    Fg = [0, 0, - G3 * G6]

    # Do the actual work here, both methods
    # X_i(t+1) = X_i(t) + SPEED * ((qMi + qli)^-1 *
    #                        (qMiX2 + qliX2 + PX2X3 + Fg) - X_i(t))
    if METHOD == 'gradient':
        new_vertex = matplus(
            vertices[i],
            matmul(
                SPEED,
                matminus(
                    matmul(
                        inverse(matplus(qMi, diag(qli, dim=3))),
                        matplus(qMiX2, matplus(qliX2, matplus(PX2X3, Fg)))
                    ),
                    vertices[i]
                )
            )
        )
    # X_i(t+1) = X_i(t) + SPEED / (ql2i + qli) *
    #                              (qMiX2 + qliX2 + PX2X3 - (qMi + qli).X_i(t))
    elif METHOD == 'seidel' or METHOD == 'fixed-point':
        new_vertex = matplus(
            vertices[i],
            matmul(
                SPEED / (ql2i + qli),
                matminus(
                    matplus(qMiX2, matplus(qliX2, matplus(PX2X3, Fg))),
                    matmul(
                        matplus(qMi, diag(qli, dim=3)),
                        vertices[i]
                    )
                )
            )
        )

    # Update displacement criterion, if needed.
    # max(squared disp) is not a really good criterion, to be improved..
    if (dotproduct(matminus(new_vertex, vertices[i]),
                   matminus(new_vertex, vertices[i]))
            > var):
        var = dotproduct(matminus(new_vertex, vertices[i]),
                         matminus(new_vertex, vertices[i]))

    if SAVE_RESULTS:
        # Branch the lists of lists to a proper GH output
        for M in matminus(qMiX2, matmul(qMi, vertices[i])):
            out['M'].Add(M / (qli + ql2i), GH_Path(iter_qs, iter, i))
        for C in matminus(qliX2, matmul(qli, vertices[i])):
            out['C'].Add(C / (qli + ql2i), GH_Path(iter_qs, iter, i))
        for p in PX2X3:
            out['P'].Add(p / (qli + ql2i), GH_Path(iter_qs, iter, i))
        for g in Fg:
            out['G'].Add(g / (qli + ql2i), GH_Path(iter_qs, iter, i))

    return var, new_vertex, out


def iterate_one_step(vertices, vertices_faces_nodes, vertices_faces,  # noqa
                     naked, fixed, qs, ql, n_cable, P6, G6, iter_qs, iter,
                     out):
    """Updates all nodes on the mesh once

    Arguments:
      vertices = list of mesh vertices, from rs.MeshVertices(mesh)
      vertices_faces_nodes = list of lists of list of nodes consituting
            the faces connected to a vertex, from orient_mesh_faces(mesh)
      vertices_faces = list of list of faces adjacent to a node, from
            rs.meshVertexFaces(mesh, i)
      naked = list of booleans, True if vertex is on a naked edge
      fixed = list of booleans, True if vertex is fixed
      qs = list of surface stress density coefficients for each face
      ql = list of cable force density, for each cable segment
      n_cable = list of list of cables segments connected to a vertex
      P6 = pressure / 6
      G6 = surfacic weight density / 6
      iter_qs = current iteration in the qs loop
      iter = current iteration number
      out = dictionnary of lists saving the forces acting on each
            vertex (qMiX2, qliX2 and PX2X3)
    Returns:
      iter = updated iteration number
      var = maximum squared displacement for this iteration
      vertices = updated positions of vertices
      out = results updated to current iteration
    """

    if METHOD == 'seidel':
        # Update the mesh in-place as we loop over each vertex
        new_vertices = vertices
    else:
        # Copy vertices to a new list, so that we do not overwrite it
        new_vertices = copy.deepcopy(vertices)
    iter += 1
    var = 0
    for i in xrange(len(vertices)):
        if not fixed[i]:
            var, new_vertices[i], out = iterate_vertex(
                i, vertices, vertices_faces_nodes, vertices_faces[i],
                naked[i], qs, ql, n_cable, P6, G6, var, iter_qs, iter, out
            )
    if METHOD == 'seidel':
        vertices == new_vertices
    else:
        vertices = copy.deepcopy(new_vertices)

    return iter, var, vertices, out


def iterate_fixed_qs(vertices, vertices_faces_nodes, vertices_faces,  # noqa
                     connec, naked, fixed, qs, ql, n_cable, P, G, iter_qs, out,
                     meshi=None):
    """Iterates the problem with qs fixed, until a pseudo-minimal
    surface is found.
    Arguments:
      vertices = list of mesh vertices, from rs.MeshVertices(mesh)
      vertices_faces_nodes = list of lists of list of nodes consituting
            the faces connected to a vertex, from orient_mesh_faces(mesh)
      vertices_faces = list of list of faces adjacent to a node, from
            rs.meshVertexFaces(mesh, i)
      connec = connectivity matrix in the mesh, from
               rs.MeshFaceVertices(mesh)
      naked = list of booleans, True if vertex is on a naked edge
      fixed = list of booleans, True if vertex is fixed
      qs = list of surface stress density coefficients for each face
      ql = list of cable force density, for each cable segment
      n_cable = list of list of cable segments attached to each vertex
      P = pressure
      G = weight density per surface area (i.e. volumic density * thickness)
      iter_qs = current iteration in the qs loop
      out = dictionnary of lists saving the forces acting on each
            vertex (qMiX2, qliX2 and PX2X3)
      meshi = save-state for graphical display
    Returns:
      vars = list of maximum squared displacement for each iteration
      vertices = updated positions of vertices
      out = results updated to current iteration
    """

    # Initialize loop
    iter = 0
    var = 2 * MAX_DISP + 1
    vars = []

    # Loop while we can, get some display if wanted
    while (iter < MAX_ITER) & (var > MAX_DISP):
        iter, var, vertices, out = iterate_one_step(
            vertices, vertices_faces_nodes, vertices_faces, naked, fixed, qs,
            ql, n_cable, P / 6, G / 6, iter_qs, iter, out
        )
        if GRAPHIC:
            rs.HideObject(meshi)
            meshi = rs.AddMesh(vertices, connec)
            print (u"itération {},"
                   u"\nDéplacement^2 maximum depuis "
                   u"l'itération précedente : {} mm²").format(iter, var)
        else:
            print var
        vars.append(var)

    return vars, vertices, out


def minimize_mesh(mesh, cables=None, fixed=None, qs=None, q_cables=None,  # noqa
                  reference=None, ql=None, n_cable=None, P=0, G=0):
    """Iterates a mesh until it is close to a minimal surface.

    Arguments:
      mesh = mesh to calculate, Rhino GUID
      cables = polylines representing the cables
      fixed = list of booleans, True if vertex is fixed
      qs = list of surface stress density coefficients for each face
      q_cables = list of force density coefficients for each cable
      reference = reference mesh for comparisons, unuseds
      ql = list of cable force density, for each cable segment
      n_cable = list of list of cable segments connected to a vertex
      P = pressure
      G = weight density per surface area (i.e. volumic density * thickness)
    Returns:
      vertices = vertices at new position
      out = dictionnary of lists saving the forces acting on each
            vertex ('M' = qMiX2, 'C' = qliX2 and 'P' = PX2X3),
            also 'S' = faces stress and 'convergence' = convergence
            criterion for each iteration.
    """

    # Initialize
    # qs defaults to 1 everywhere
    # meshi remembers current mesh if we need to hide it in Rhino
    # vertices strips RhinoCommon's 3D points to a bare 3-vector
    # reference defaults to the initial mesh
    # vertices_faces is taken from rs.MeshVertexFaces(mesh, i)
    # vertices_faces_nodes is reordered from Rhino's mesh
    # naked is extracted from Rhino
    # fixed defaults to naked if edges are fixed, otherwise defined in
    #                   mmh.define_cables
    # ql, n_cable are defined by mmh.define_cables if necessary
    # q_cables defaults to 1
    # out is just empty GH trees at first

    if not qs:
        qs = [1 for i in range(rs.MeshFaceCount(mesh))]
    elif type(qs) != list:
        qs = [qs for i in range(rs.MeshFaceCount(mesh))]
        print qs
    if GRAPHIC or SHOW_RESULT:
        meshi = mesh
    old_vertices = rs.MeshVertices(mesh)
    vertices = [[old_vertices[i][0], old_vertices[i][1], old_vertices[i][2]]
                for i in range(len(old_vertices))]
    if not reference:
        reference = vertices
    vertices_faces = [rs.MeshVertexFaces(mesh, i)
                      for i in range(rs.MeshVertexCount(mesh))]
    vertices_faces_nodes = mmh.orient_mesh_faces(mesh)
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
    out = {'M': DataTree[object](),
           'C': DataTree[object](),
           'P': DataTree[object](),
           'G': DataTree[object](),
           'S': [0 for i in range(rs.MeshFaceCount(mesh))],
           'convergence': DataTree[object](),
           }

    # Initialize loop
    iter_qs = 0
    dev_sigma = 2 * MAX_DEV_SIGMA + 1

    # Loop while qs needs to be updated
    while (iter_qs < MAX_ITER_QS) & (dev_sigma > MAX_DEV_SIGMA):
        vars, vertices, out = iterate_fixed_qs(
            vertices, vertices_faces_nodes, vertices_faces, connec, naked,
            fixed, qs, ql, n_cable, P, G, iter_qs, out
        )
        # Branch the list of lists to a proper GH output
        for i, var in enumerate(vars):
            out['convergence'].Add(var, GH_Path(*[0, iter_qs]))
        iter_qs += 1
        if iter_qs < MAX_ITER_QS:
            dev_sigma, qs, out['S'] = mmh.update_qs(mesh, qs)
        else:
            out['S'] = mmh.update_qs(mesh, copy.deepcopy(qs))[2]
        if DEBUG:
            print qs
        print '-' * 20

    if SHOW_RESULT:
        rs.HideObject(meshi)
        rs.AddMesh(vertices, connec)

    return vertices, out
