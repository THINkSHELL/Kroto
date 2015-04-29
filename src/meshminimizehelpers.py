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

"""Defines some helper functions for the meshminimize module
"""

import rhinoscriptsyntax as rs
import vectorworks as vw
import meshminimize as mm


def redraw_off(func, *args, **kargs):
    """Decorator function turning off Redraw in Rhino.

    Parameters:
      func = function to be decorated
      *args, **kargs = arguments
    Returns:
      wrapper = decorated function
    """
    def wrapper(*args, **kargs):
        rs.EnableRedraw(False)
        result = func(*args, **kargs)
        rs.EnableRedraw(True)
        rs.Redraw()
        return result
    return wrapper


@redraw_off
def upward_face(n, x1, x2, x3):
    """Returns True if the face orientation defined by the order of the
    passed in points (screw rule) is the same as that of the mesh,
    passed by the normal n.

    Parameters:
        n = Vector normal of the face
            n would normally be taken from rs.MeshFaceNormals(mesh)
            which gives a consistent orientation to the mesh over all
            faces, if normals have unified.
        x1, x2, x3 = 3D Points representing the triangular face
    Returns:
        True or False
    """
    v2 = rs.VectorCreate(x2, x1)
    v3 = rs.VectorCreate(x3, x1)

    # Using rs's method for vector computation is slower than vw's, but
    # it takes into account nasty tolerance things.
    vv = rs.VectorCrossProduct(v2, v3)
    a = rs.VectorDotProduct(vv, n)
    return a > 0


@redraw_off  # noqa
def define_cables(cables, q_cables, vertices, naked, fixed):
    """Defines the cables list strucure from the Rhino geometry.
    Connects the mesh vertices lying on a polyline together. The
    polyline vertices are not considered as ends for the
    mm.FIXED_CABLE_ENDS option.

    Arguments:
      cables = list of Rhino polylines representing the cables
      q_cables = list of force density coefficients for each cable
      vertices = list of mesh vertices
      naked = list of naked mesh vertices, from rs.MeshNakedVertice(mesh)
      fixed = list of fixed vertices
    Returns:
      ql = list of list of force density coef for each cable segment
           connected to the vertices
      n_cable = number of cable segments connected to each vertices
      fixed = updated list of fixed vertices
    """

    # Make sure we know what "point on cable" means
    tol = rs.UnitAbsoluteTolerance()

    # Initialize
    # vcab = list of list of vertices on each cable
    # ql = list of list of the connected cables force densities,
    #      for each vertex
    # n_cable = connectivity matrix of the cables
    #         = list of list of connected vertices to each vertex
    vcab = [[] for i in cables]
    ql = [[] for i in vertices]
    n_cable = [[] for i in vertices]

    for v, vertex in enumerate(vertices):
        # Only consider naked edge vertices, if a cable is in the middle
        #  of a mesh then these edges have to be split
        if naked[v]:
            for i, cable in enumerate(cables):
                v_proj_cab = rs.EvaluateCurve(
                    cable,
                    rs.CurveClosestPoint(cable, vertex)
                )
                dis_v_to_cab = rs.Distance(v_proj_cab, vertex)

                if dis_v_to_cab < tol:
                    # Vertex is on cable i, save it to vcab[i]
                    vcab[i].append(v)
                    dis_v_to_ends = min(
                        rs.Distance(rs.CurveEndPoint(cable), vertex),
                        rs.Distance(rs.CurveStartPoint(cable), vertex)
                    )
                    if mm.FIXED_CABLE_ENDS and dis_v_to_ends < tol:
                        # Vertex is on a cable end, fix it if needed
                        if mm.DEBUG:
                            rs.AddPoint(vertex)
                        fixed[v] = True

    for i, cable in enumerate(cables):
        if vcab[i] == []:
            raise ValueError(
                'Cable #%i appears to be too far off the mesh' % i
            )
        # Sort vertices along cable, successive vertices will be linked
        vcab[i].sort(key=lambda v: rs.CurveClosestPoint(cable, vertices[v]))
        # If cable is closed, re-add first vertex at the end
        if rs.Distance(rs.CurveStartPoint(cable),
                       rs.CurveEndPoint(cable)) < tol:
            vcab[i].append(vcab[i][0])

        # Construct connectivity matrix for current cable
        for j in range(len(vcab[i])):
            if not fixed[vcab[i][j]]:
                if j:
                    n_cable[vcab[i][j]].append(vcab[i][j - 1])
                    ql[vcab[i][j]].append(q_cables[i])
                if j - len(vcab[i]) + 1:
                    n_cable[vcab[i][j]].append(vcab[i][j + 1])
                    ql[vcab[i][j]].append(q_cables[i])

    return ql, n_cable, fixed


@redraw_off
def orient_mesh_faces(mesh):
    """Orients faces around the nodes in a mesh to a consistant order
    and normal direction. Roughly equivalent to rs.MeshFaceVertices, but
    we control the list order.

    Arguments:
      mesh = the mesh in RhinoCommon type
    Returns:
      vertex_faces list of faces adjacent to a node
    """

    rs.EnableRedraw(False)
    normals = rs.MeshFaceNormals(mesh)
    vertices = rs.MeshVertices(mesh)
    connec = rs.MeshFaceVertices(mesh)
    vertex_faces = {}
    for i, vertex in enumerate(vertices):
        for j, face in enumerate(rs.MeshVertexFaces(mesh, i)):
            # Find the three distinct points, 4th is redundant in
            # triangular meshes and remove current point from the set.
            # Might break the original ordering, but we will take care
            # of that our own way afterwards.
            others = list(set(connec[face]) - {i, })
            [x2, x3] = [vertices[n] for n in others]
            if upward_face(normals[face], vertex, x2, x3):
                vertex_faces[(i, j)] = [i, others[0], others[1]]
                if mm.DEBUG:
                    print 'not flip'
            else:
                vertex_faces[(i, j)] = [i, others[1], others[0]]
                if mm.DEBUG:
                    print 'flip'
    rs.EnableRedraw
    return vertex_faces


@redraw_off
def mesh_distance(vertices, objective):
    """Evauluates the distance between two set of vertices representing
    the same mesh in different positions.  The distance is just a max of
    the squared length between the two vertices at corresponding indices
    in the list, so the meshes should look similar from the beginning if
    we want this to mean something.

    Arguments:
      vertices, objectives = list of vertices, from rs.MeshVertices(mesh)
    Returns:
      max squared distance between two vertices at the same index
    """

    distance = 0
    for i in range(len(vertices)):
        temp = vw.dist3(objective[i], vertices[i])
        if temp > distance:
            distance = temp
    return distance


@redraw_off
def mesh_closest_vertices(mesh, points):
    """Finds the vertices of a mesh that are within the document's tolerance
    of one of the points in the list.
    Arguments:
        mesh = a Rhino mesh
        points = a list of Rhino points
    Returns:
        fixed = a list of booleans, True if the vertex is close to one of
                the points
    """

    import rhinoscriptsyntax as rs

    tol = rs.UnitAbsoluteTolerance()
    vertices = rs.MeshVertices(mesh)
    fixed = [False for i in vertices]

    for i, vertex in enumerate(vertices):
        for point in points:
            if rs.Distance(point, vertex) < tol:
                fixed[i] = True
                break
    return fixed


def update_qs(mesh, qs):
    """Updates the surface stress density coefficients to reach a more
    uniform stress over the surface. Computes surface stresses at the
    same time.
    Arguments:
      mesh = the Rhino mesh
      qs = the list current values of qs for each face
    Returns:
      dev_sigma = maximum deviation to mean-value of the surface stress
      qs = updated qs list
      sigma = surface stresses
    """

    faces = rs.MeshFaceVertices(mesh)
    vertices = rs.MeshVertices(mesh)
    mean_sigma = 0
    dev_sigma = 0
    sigma = [0 for i in faces]
    n_faces = len(sigma)

    # Compute mean stress value
    for i, face in enumerate(faces):
        x1 = vertices[face[0]]  # = X_1j, 3D point, [3x1] vector
        x2 = vertices[face[1]]  # = X_2j, 3D point, [3x1] vector
        x3 = vertices[face[2]]  # = X_3j, 3D point, [3x1] vector
        x12 = vw.vecminus3(x2, x1)
        x23 = vw.vecminus3(x3, x1)
        face_area = .5 * vw.norm3(vw.crossproduct3(x12, x23))
        sigma[i] = qs[i] * face_area
        mean_sigma += sigma[i] / n_faces

    # Find the maximum deviation to the mean value
    dev_sigma = max([abs(sig - mean_sigma) for sig in sigma])

    # If above tolerance, update
    if dev_sigma > mm.MAX_DEV_SIGMA:
        for i, sigma_i in enumerate(sigma):
            qs[i] = qs[i] * mean_sigma / sigma_i

    return dev_sigma, qs, sigma
