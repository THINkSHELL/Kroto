"""Defines some helper functions for the meshminimize module
"""

import rhinoscriptsyntax as rs
import vectorworks as vw
import meshminimize as mm


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


def define_cables(cables, q_cables, vertices, naked, fixed):  # noqa
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
    # temp = list of list of vertices on each cable
    # ql = list of list of the connected cables force densities,
    #      for each vertex
    # n_cable = connectivity matrix of the cables
    #         = list of list of connected vertices to each vertex
    temp = [[] for i in cables]
    ql = [[] for i in vertices]
    n_cable = [[] for i in vertices]

    # Only consider naked edge vertices, if a cable is in the middle of
    # a mesh then these edges have to be split
    for v, vertex in enumerate(vertices):

        if naked[v]:
            for i, cable in enumerate(cables):

                # Vertex is on cable i, save it to temp[i]
                if (rs.Distance(
                        rs.EvaluateCurve(
                            cable,
                            rs.CurveClosestPoint(cable, vertex)),
                        vertex) < tol):
                    temp[i].append(v)

                    # Vertex is on a cable end, fix it if needed
                    if (mm.FIXED_CABLE_ENDS and
                            min(rs.Distance(rs.CurveEndPoint(cable), vertex),
                                rs.Distance(rs.CurveStartPoint(cable), vertex))
                            < tol):
                        if mm.DEBUG:
                            rs.AddPoint(vertex)
                        fixed[v] = True

    for i, cable in enumerate(cables):

        if temp[i] == []:
            raise ValueError('Cable #%i is to be too far off the mesh' % i)

        # Sort vertices along cable, successive vertices will be linked
        temp[i].sort(key=lambda v: rs.CurveClosestPoint(cable, vertices[v]))

        # If cable is closed, re-add first vertex at the end
        if rs.Distance(rs.CurveStartPoint(cable),
                       rs.CurveEndPoint(cable)) < tol:
            print rs.Distance(rs.CurveStartPoint(cable),
                              rs.CurveEndPoint(cable))
            temp[i].append(temp[i][0])

        for j in range(len(temp[i])):
            if not fixed[temp[i][j]]:
                if j:
                    n_cable[temp[i][j]].append(temp[i][j - 1])
                    ql[temp[i][j]].append(q_cables[i])
                if (j - len(temp[i]) + 1):
                    n_cable[temp[i][j]].append(temp[i][j + 1])
                    ql[temp[i][j]].append(q_cables[i])

    return ql, n_cable, fixed


def orient_mesh_faces(mesh):
    """Orients faces around the nodes in a mesh to a consistant order
    and normal direction. Roughly equivalent to rs.MeshFaceVertices, but
    we control the list order.

    Arguments:
      mesh = the mesh in RhinoCommon type
    Returns:
      vertex_faces list of faces adjacent to a node
    """

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
    return vertex_faces


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
        if (vw.dotproduct(vw.matminus(objective[i], vertices[i]),
                          vw.matminus(objective[i], vertices[i]))
                > distance):
            distance = vw.dotproduct(vw.matminus(objective[i], vertices[i]),
                                     vw.matminus(objective[i], vertices[i]))
    return distance


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
        x12 = vw.matminus(x2, x1)
        x23 = vw.matminus(x3, x1)
        n = vw.crossproduct(x12, x23)
        face_area = .5 * vw.dotproduct(n, n) ** .5
        sigma[i] = qs[i] * face_area
        mean_sigma += sigma[i] / n_faces

    # Find the maximum deviation to the mean value
    dev_sigma = max([abs(sig - mean_sigma) for sig in sigma])

    # If above tolerance, update
    if dev_sigma > mm.MAX_DEV_SIGMA:
        for i, sigma_i in enumerate(sigma):
            qs[i] = qs[i] * mean_sigma / sigma_i

    return dev_sigma, qs, sigma
