#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Kroto Membrane Form-finding.
#  Copyright © 2015, Thinkshell & Laboratoire Navier, ENPC.
#
#  This file is part of Kroto.
#
#  Kroto is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  Kroto is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with Kroto. If not, see http://www.gnu.org/licenses/.
#

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
from imp import reload
reload(vw)

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
DEBUG = False
GRAPHIC = False
SHOW_RESULT = False
SAVE_RESULTS = False
MAX_DISP = 0.01
MAX_ITER = 10
MAX_ITER_QS = 2
MAX_DEV_SIGMA = 1
SPEED = 1
METHOD = 'seidel'
FIXED_CABLE_ENDS = True


def iterate_vertex(i, vertices, vertices_faces_nodes, vertex_faces, naked, qs,
                   ql, n_cable, p6, g6, res, iter_qs, iter, out):
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
      p6 = pressure / 6
      g6 = surfacic weight density / 6
      res = the current maximum displacement in the iteration
      iter_qs = current iteration in the qs loop
      iter = current iteration in the stresses loop
      out = dictionnary of lists saving the forces acting on each
            vertex (q_mix2, q_lix2 and p_x2x3)
    Returns
      res = the updated maximum displacement
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
    # x_2j, x_3j = points 2 and 3 of the face #j around the vertex i,
    #              [3x1] vector
    # x_2j3j = x_3j - x_2j = vector from 2j to 3j, [3x1] vector

    # M_(i,j) = (x_2j3j . x_2j3j)*Id - (x_2j3j x x_2j3j)
    #           [3,3] matrix representing the dependant faces areas around
    #           the vertex i

    # qs_j  = face number j surface stress density coefficient
    # ql_j  = cable segment number j force density coefficient (for points
    #         in the middle of a cable, count each side once)
    # p6    = pressure / 6 (uniform scalar at the moment)
    # g6    = surfacic weight density / 6 (uniform scalar at the moment)

    # q_mi   = sum(j = [1, m_i]; qs_j * M_(i,j))
    #       = local stiffness matrix, [3x3] matrix
    # q_mix2 = sum(j = [1, m_i]; qs_j * M_(i,j) * x_2j)
    #       = local membrane forces on the vertex, [3x1] vector
    # ql2i  = sum(j = [1, m_i]; qs_j * l_ij**2)
    #       = local membrane stifness coefficient, scalar
    # qli   = sum(j = [1, n_i]; ql_j)
    #       = local cable stifness coefficient, scalar
    # q_lix2 = sum(j = [1, n_i]; ql_j * M_(i,j) * x_2j)
    #       = local cable forces on the vertex, [3x1] vector
    # p_x2x3 = sum(j = [1, m_i]; p/6 * (x_2j/\x_3j + x_2j3j/\x_i))
    #       = pressure membrane forces around the vertex, [3x1] vector
    # g3    = sum(j = [1, m_i]; norm(x_12j/\x_13j))
    #       = u3 term / g6 of the gravity forces around the vertex, scalar
    # f_g    = [0, 0, - g3 * g6] = gravity forces, [3x1] vector

    # First we initialize the intermediate vectors for calculation
    q_mi = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    q_mix2 = [0, 0, 0]
    ql2i = 0
    if ql:
        qli = sum(ql[i])
    else:
        qli = 0
    q_lix2 = [0, 0, 0]
    p_x2x3 = [0, 0, 0]
    g3 = 0

    j = 0

    if DEBUG:
        print '####'
    if DEBUG:
        print vertices[i]

    # Iterate for faces j around current vertex i
    while vertices_faces_nodes.get((i, j), 0):  # Face is in vertex_faces table
        x2 = vertices[vertices_faces_nodes[(i, j)][1]]  # = x_2j,3D point,[3x1]
        x3 = vertices[vertices_faces_nodes[(i, j)][2]]  # = x_3j,3D point,[3x1]
        x23 = vecminus3(x3, x2)  # = vector x_2j3j, [3x1] vector
        if g6:
            x12 = vecminus3(x2, vertices[i])
            x13 = vecminus3(x3, vertices[i])
        qij = qs[vertex_faces[j]]  # = qs_j = Surface stress density coef

        # q_mij = qs_j * M_(i,j)
        # M_(i,j) = (x_2j3j . x_2j3j)*Id - (x_2j3j x x_2j3j)
        q_mij = scalmatmul3(
            qij,
            matminus3(
                diag3(dotproduct3(x23, x23)),
                veckronproduct3(x23, x23)
            )
        )

        # update q_mi and q_mix2 running sums
        q_mi = matplus3(q_mi, q_mij)
        q_mix2 = vecplus3(q_mix2, matvecmul3(q_mij, x2))

        # ql2i = sum(j = [1, m_i]; qs_j * l_ij**2), running sum update, scalar
        ql2i += qij * dotproduct3(x23, x23)

        # p_x2x3 = sum(j = [1, m_i]; p/6 * (x_2j/\x_3j + x_2j3j/\x_i),
        # running sum update, [3,1] vector
        if p6:
            p_x2x3 = vecplus3(
                p_x2x3,
                scalvecmul3(p6, crossproduct3(x2, x3))
            )
            # The contour term is non-zero only for edge nodes
            if naked:
                p_x2x3 = vecplus3(
                    p_x2x3,
                    scalvecmul3(p6, crossproduct3(x23, vertices[i]))
                )

        # g3 = sum(j = [1, m_i]; norm(x_12j/\x_13j)),
        # running sum update, scalar
        if g6:
            g3 += norm3(crossproduct3(x12, x13))

        if DEBUG:
            print '\n'.join([
                str((i, j)), str(x2), str(x3), str(x23), str(ql),
                str(q_mij), matvecmul3(q_mij, x2), str(ql2i), '---'
            ])
        j += 1
        # END while

    if DEBUG:
        print '\n'.join([q_mi, q_mix2, p_x2x3])

    # Iterate for cables segments j around vertex i
    if n_cable:
        for j, v in enumerate(n_cable[i]):
            # q_lix2 = q_lix2 = sum( j = [1, n_i]; ql_j * M_(i,j) * x_2j ),
            # Running sum update, [3x1] vector.
            q_lix2 = vecplus3(q_lix2, scalvecmul3(ql[i][j], vertices[v]))

    # Build the gravity forces vector
    # f_g = [0, 0, - g3 * g6]
    f_g = [0, 0, - g3 * g6]

    # Do the actual work here, both methods
    # x_i(t+1) = x_i(t) + SPEED * ((q_mi + qli)**-1 *
    #                        (q_mix2 + q_lix2 + p_x2x3 + f_g) - x_i(t))
    if METHOD == 'gradient':
        new_vertex = vecplus3(
            vertices[i],
            scalvecmul3(
                SPEED,
                vecminus3(
                    matvecmul3(
                        inverse3(matplus3(q_mi, diag3(qli))),
                        vecplus34(q_mix2, q_lix2, p_x2x3, f_g)
                    ),
                    vertices[i]
                )
            )
        )
    # x_i(t+1) = x_i(t)
    #            + SPEED / (ql2i + qli)
    #              * (q_mix2 + q_lix2 + p_x2x3 - (q_mi + qli).x_i(t))
    elif METHOD == 'seidel' or METHOD == 'fixed-point':
        new_vertex = vecplus3(
            vertices[i],
            scalvecmul3(
                SPEED / (ql2i + qli),
                vecminus3(
                    vecplus34(q_mix2, q_lix2, p_x2x3, f_g),
                    matvecmul3(
                        matplus3(q_mi, diag3(qli)),
                        vertices[i]
                    )
                )
            )
        )

    # Update displacement criterion, if needed.
    temp = dist3(new_vertex, vertices[i])
    if temp > res:
        res = temp

    if SAVE_RESULTS:
        # Branch the lists of lists to a proper GH output
        for M in vecminus3(q_mix2, matvecmul3(q_mi, vertices[i])):
            out['M'].Add(M / (qli + ql2i), GH_Path(iter_qs, iter, i))
        for C in vecminus3(q_lix2, scalvecmul3(qli, vertices[i])):
            out['C'].Add(C / (qli + ql2i), GH_Path(iter_qs, iter, i))
        for p in p_x2x3:
            out['P'].Add(p / (qli + ql2i), GH_Path(iter_qs, iter, i))
        for g in f_g:
            out['G'].Add(g / (qli + ql2i), GH_Path(iter_qs, iter, i))

    return res, new_vertex, out


def iterate_one_step(vertices, vertices_faces_nodes, vertices_faces, naked,
                     fixed, qs, ql, n_cable, p6, g6, iter_qs, iter, out):
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
      p6 = pressure / 6
      g6 = surfacic weight density / 6
      iter_qs = current iteration in the qs loop
      iter = current iteration number
      out = dictionnary of lists saving the forces acting on each
            vertex (q_mix2, q_lix2 and p_x2x3)
    Returns:
      iter = updated iteration number
      res = maximum squared displacement for this iteration
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
    res = 0

    for i in xrange(len(vertices)):
        if not fixed[i]:
            res, new_vertices[i], out = iterate_vertex(
                i, vertices, vertices_faces_nodes, vertices_faces[i],
                naked[i], qs, ql, n_cable, p6, g6, res, iter_qs, iter, out
            )
    if METHOD == 'seidel':
        vertices == new_vertices
    else:
        vertices = copy.deepcopy(new_vertices)

    return iter, res, vertices, out


def iterate_fixed_qs(vertices, vertices_faces_nodes, vertices_faces, connec,
                     naked, fixed, qs, ql, n_cable, p, g, iter_qs, out,
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
      p = pressure
      g = weight density per surface area (i.e. volumic density * thickness)
      iter_qs = current iteration in the qs loop
      out = dictionnary of lists saving the forces acting on each
            vertex (q_mix2, q_lix2 and p_x2x3)
      meshi = save-state for graphical display
    Returns:
      res_list = list of maximum squared displacement for each iteration
      vertices = updated positions of vertices
      out = results updated to current iteration
    """

    # Initialize loop
    iter = 0
    res = 2 * MAX_DISP + 1
    res_list = []

    # Loop while we can, get some display if wanted
    while (iter < MAX_ITER) & (res > MAX_DISP):
        iter, res, vertices, out = iterate_one_step(
            vertices, vertices_faces_nodes, vertices_faces, naked, fixed, qs,
            ql, n_cable, p / 6, g / 6, iter_qs, iter, out
        )
        res_list.append(res)

        if GRAPHIC:
            rs.HideObject(meshi)
            meshi = rs.AddMesh(vertices, connec)
            print (u"Iteration number {},"
                   u"\nMaximum displacement since "
                   u"previous iteration : {} mm").format(iter, res)
        else:
            print res
        if SAVE_RESULTS:
            out['A'].Add(rs.MeshArea(rs.AddMesh(vertices, connec))[1],
                         GH_Path(iter_qs))

    return res_list, vertices, out


def minimize_mesh(mesh, cables=None, fixed=None, qs=None, q_cables=None,
                  reference=None, ql=None, n_cable=None, p=0, g=0):
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
      p = pressure
      g = weight density per surface area (i.e. volumic density * thickness)
    Returns:
      vertices = vertices at new position
      out = dictionnary of lists saving the forces acting on each
            vertex ('M' = q_mix2, 'C' = q_lix2 and 'p' = p_x2x3),
            also 'S' = faces stress and 'convergence' = convergence
            criterion for each iteration.
    """

    # INITIALIZE
    # qs defaults to 1 everywhere
    if not qs:
        qs = [1 for i in range(rs.MeshFaceCount(mesh))]
    elif type(qs) != list:
        qs = [qs for i in range(rs.MeshFaceCount(mesh))]
        print qs

    # meshi remembers current mesh if we need to hide it in Rhino
    if GRAPHIC or SHOW_RESULT:
        meshi = mesh
    old_vertices = rs.MeshVertices(mesh)

    # vertices strips RhinoCommon's 3D points to a bare 3-vector
    vertices = [[old_vertices[i][0], old_vertices[i][1], old_vertices[i][2]]
                for i in range(len(old_vertices))]

    # reference defaults to the initial mesh
    if not reference:
        reference = vertices

    # vertices_faces is taken from rs.MeshVertexFaces(mesh, i)
    vertices_faces = [rs.MeshVertexFaces(mesh, i)
                      for i in range(rs.MeshVertexCount(mesh))]

    # vertices_faces_nodes is reordered from Rhino's mesh
    vertices_faces_nodes = mmh.orient_mesh_faces(mesh)

    # naked is extracted from Rhino
    naked = rs.MeshNakedEdgePoints(mesh)
    connec = rs.MeshFaceVertices(mesh)

    # fixed defaults to naked if edges are fixed, otherwise defined in
    #                   mmh.define_cables
    if not fixed:
        if not cables:
            fixed = naked
        else:
            fixed = [False for i in vertices]

    # ql, n_cable are defined by mmh.define_cables if necessary
    ql = None
    n_cable = None

    # q_cables defaults to 1
    if cables:
        if not q_cables:
            q_cables = [1 for i in cables]
        if not (ql and n_cable and fixed):
            ql, n_cable, fixed = mmh.define_cables(cables, q_cables,
                                                   old_vertices, naked, fixed)

    # out is just empty GH trees at first
    out = {'M': DataTree[object](),
           'C': DataTree[object](),
           'P': DataTree[object](),
           'G': DataTree[object](),
           'A': DataTree[object](),
           'S': [0 for i in range(rs.MeshFaceCount(mesh))],
           'convergence': DataTree[object](),
           }

    # Initialize loop
    iter_qs = 0
    dev_sigma = 2 * MAX_DEV_SIGMA + 1

    # Loop while qs needs to be updated
    while (iter_qs < MAX_ITER_QS) & (dev_sigma > MAX_DEV_SIGMA):
        res_list, vertices, out = iterate_fixed_qs(
            vertices, vertices_faces_nodes, vertices_faces, connec, naked,
            fixed, qs, ql, n_cable, p, g, iter_qs, out
        )
        # Branch the list of lists to a proper GH output
        for i, res in enumerate(res_list):
            out['convergence'].Add(res, GH_Path(*[0, iter_qs]))
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
