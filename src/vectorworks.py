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

"""Provides simple vector and matrix calculus functions.
"""


def dotproduct(x1, x2):
    """Computes the cross-product of 2 vectors.
    """
    return sum([x1[i] * x2[i] for i in range(len(x1))])


def dotproduct3(x1, x2):
    """Faster cross-product of 2 vectors, in dimension 3 ONLY.
    """
    return x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2]


def crossproduct3(x1, x2):
    """Computes the cross-product of 2 vectors, in dimension 3 ONLY.
    """
    return [
        x1[1] * x2[2] - x1[2] * x2[1],
        x1[2] * x2[0] - x1[0] * x2[2],
        x1[0] * x2[1] - x1[1] * x2[0]
    ]


def inverse3(a):
    """Computes the inverse of a matrix, in dimension 3 ONLY.
    Returns None if the determinant is 0.
    """
    det = (
        a[0][0] * a[1][1] * a[2][2] +
        a[1][0] * a[2][1] * a[0][2] +
        a[2][0] * a[0][1] * a[1][2] -
        a[0][2] * a[1][1] * a[2][0] -
        a[1][2] * a[2][1] * a[0][0] -
        a[2][2] * a[0][1] * a[1][0]
    )
    if det:
        res = [
            [
                a[1][1] * a[2][2] / det - a[1][2] * a[2][1] / det,
                a[0][2] * a[2][1] / det - a[0][1] * a[2][2] / det,
                a[0][1] * a[1][2] / det - a[0][2] * a[1][1] / det
            ],
            [
                a[1][2] * a[2][0] / det - a[1][0] * a[2][2] / det,
                a[0][0] * a[2][2] / det - a[0][2] * a[2][0] / det,
                a[0][2] * a[1][0] / det - a[0][0] * a[1][2] / det
            ],
            [
                a[1][0] * a[2][1] / det - a[1][1] * a[2][0] / det,
                a[0][1] * a[2][0] / det - a[0][0] * a[2][1] / det,
                a[0][0] * a[1][1] / det - a[0][1] * a[1][0] / det
            ]
        ]
    else:
        res = None
    return res


def matplus(a, b):
    """Computes the term-by-term sum of 2 array-like objects.
    Accepts any combination of scalars, vectors and matrices,
    as long as their dimensions match.
    """
    if type(a) != list and type(b) != list:
        res = a + b
    elif type(a[0]) != list and type(b[0]) != list:
        res = [a[i] + b[i] for i in range(len(a))]
    else:
        try:
            res = [[a[i][j] + b[i][j] for j in range(len(b[0]))]
                   for i in range(len(a))]
        except TypeError:
            raise TypeError(
                "Matrices dimensions not compatible: %s and %s" % (a, b))
    return res


def matplus3(a, b):
    """Faster term-by-term matrix sum, in dimension 3 ONLY.
    """
    return [
        [a[0][0] + b[0][0], a[0][1] + b[0][1], a[0][2] + b[0][2]],
        [a[1][0] + b[1][0], a[1][1] + b[1][1], a[1][2] + b[1][2]],
        [a[2][0] + b[2][0], a[2][1] + b[2][1], a[2][2] + b[2][2]]
    ]


def vecplus3(x1, x2):
    """Faster term-by-term vector sum, in dimension 3 ONLY.
    """
    return [x1[0] + x2[0], x1[1] + x2[1], x1[2] + x2[2]]


def vecplus34(x1, x2, x3, x4):
    """Faster term-by-term 4 vector sum, in dimension 3 ONLY.
    """
    return [
        x1[0] + x2[0] + x3[0] + x4[0],
        x1[1] + x2[1] + x3[1] + x4[1],
        x1[2] + x2[2] + x3[2] + x4[2]
    ]


def matminus(a, b):
    """Computes the term-by-term difference of 2 array-like objects.
    Accepts any combination of scalars, vectors and matrices,
    as long as their dimensions match.
    """
    if type(a) != list and type(b) != list:
        res = a - b
    elif type(a[0]) != list and type(b[0]) != list:
        res = [a[i] - b[i] for i in range(len(a))]
    else:
        try:
            res = [[a[i][j] - b[i][j] for j in range(len(b[0]))]
                   for i in range(len(a))]
        except TypeError:
            raise TypeError(
                "Matrices dimensions not compatible: %s and %s" % (a, b))
    return res


def matminus3(a, b):
    """Faster term-by-term matrix difference, in dimension 3 ONLY.
    """
    return [
        [a[0][0] - b[0][0], a[0][1] - b[0][1], a[0][2] - b[0][2]],
        [a[1][0] - b[1][0], a[1][1] - b[1][1], a[1][2] - b[1][2]],
        [a[2][0] - b[2][0], a[2][1] - b[2][1], a[2][2] - b[2][2]]
    ]


def vecminus3(x1, x2):
    """Faster term-by-term vector difference, in dimension 3 ONLY.
    """
    return [x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2]]


def matmul(a, b):
    """Computes the classical product of array-like objects.
    Accepts any combination of scalars, vectors and matrices,
    as long as their dimensions match in the matrix product sense.
    """
    if type(a) != list and type(b) != list:
        res = a * b
    elif type(a) != list:
        if type(b[0]) == list:
            res = [[a * b[i][j] for j in range(len(b[0]))]
                   for i in range(len(b))]
        else:
            res = [a * b[i] for i in range(len(b))]
    elif type(b) != list:
        res = matmul(b, a)
    elif type(a[0]) == list and type(b[0]) == list:
        res = [[sum(a[i][k] * b[k][j] for k in range(len(b)))
               for j in range(len(b[0]))] for i in range(len(a))]
    elif type(a[0]) == list:
        res = [sum(a[i][k] * b[k] for k in range(len(b)))
               for i in range(len(a))]
    elif type(b[0]) == list:
        res = [sum(a[k] * b[k][j] for k in range(len(b)))
               for j in range(len(b[0]))]
    else:
        res = dotproduct(a, b)
    return res


def matvecmul3(a, x):
    """Faster matrix * vector multiplication, in dimension 3 ONLY.
    """
    return [
        x[0] * a[0][0] + x[1] * a[0][1] + x[2] * a[0][2],
        x[0] * a[1][0] + x[1] * a[1][1] + x[2] * a[1][2],
        x[0] * a[2][0] + x[1] * a[2][1] + x[2] * a[2][2]
    ]


def scalmatmul3(s, a):
    """Faster scalar * matrix multiplication, in dimension 3 ONLY.
    """
    return [
        [s * a[0][0], s * a[0][1], s * a[0][2]],
        [s * a[1][0], s * a[1][1], s * a[1][2]],
        [s * a[2][0], s * a[2][1], s * a[2][2]]
    ]


def scalvecmul3(s, x):
    """Faster scalar * vector multiplication, in dimension 3 ONLY.
    """
    return [s * x[0], s * x[1], s * x[2]]


def veckronproduct(x1, x2):
    """Computes the kronecker product matrix of 2 vectors.
    """
    return [[x1[i] * x2[j] for j in range(len(x2))] for i in range(len(x1))]


def veckronproduct3(x1, x2):
    """Faster vector kronecker product, in dimension 3 ONLY.
    """
    return [
        [x1[0] * x2[0], x1[0] * x2[1], x1[0] * x2[2]],
        [x1[1] * x2[0], x1[1] * x2[1], x1[1] * x2[2]],
        [x1[2] * x2[0], x1[2] * x2[1], x1[2] * x2[2]]
    ]


def transpose(a):
    """Returns the input matrix transposed.
    """
    if type(a[0]) == list:
        res = [[a[i][j] for i in range(len(a))] for j in range(len(a[0]))]
    else:
        res = [[a[i]] for i in range(len(a))]
    return res


def transpose3(a):
    """Faster matrix transposition, in dimension 3 ONLY.
    """
    return [
        [a[0][0], a[1][0], a[2][0]],
        [a[0][1], a[1][1], a[2][1]],
        [a[0][2], a[1][2], a[2][2]]
    ]


def diag(x, dim=None):
    """Returns a diagonal matrix.
    Arguments:
      x = a vector or a scalar. If x is a vector, then dim is omitted
      dim (optional) = dimension of the matrix to generate. Only
                       used if x is a scalar.
    """
    if type(x) == list:
        return [[(i == j) * x[i] for i in range(len(x))]
                for j in range(len(x))]
    elif dim:
        return [[(i == j) * x for i in range(dim)] for j in range(dim)]
    return None


def diag3(x):
    """Faster diagonal matrix from a scalar, in dimension 3 ONLY.
    """
    return [
        [x, 0, 0],
        [0, x, 0],
        [0, 0, x]
    ]


def diagvec3(x):
    """Faster diagonal matrix from a vector, in dimension 3 ONLY.
    """
    return [
        [x[0], 0, 0],
        [0, x[1], 0],
        [0, 0, x[2]]
    ]


def norm3(x):
    """Faster norm computation, in dimension 3 ONLY.
    """
    return (x[0] ** 2 + x[1] ** 2 + x[2] ** 2) ** 0.5


def dist3(x1, x2):
    """Faster distance computation, in dimension 3 ONLY.
    """
    return (
        (x1[0] - x2[0]) ** 2 +
        (x1[1] - x2[1]) ** 2 +
        (x1[2] - x2[2]) ** 2
    ) ** 0.5


if __name__ == '__main__':
    a = [[1, 2, 3], [4, 5, 6], [7, 8, 10]]
    b = [[1, 2], [3, 4]]
    x1 = [1, 2, 3]
    x2 = [4, 5, 6]
    print a
    print transpose(a)
    print dotproduct(x1, x2)
    print crossproduct3(x1, x2)
    print veckronproduct3(x1, x2)
    print inverse3(a)
    print matmul(a, x1)
    print matmul(x2, a)
    print transpose(x1)
    print matmul(3, x1)
    print matmul(a, 3)
    print diag(x1)
    print diag(5, dim=4)
