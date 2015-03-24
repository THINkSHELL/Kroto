id = [[1, 0, 0],
      [0, 1, 0],
      [0, 0, 1]]


def dotproduct(x1, x2):
    return sum([x1[i]*x2[i] for i in range(len(x1))])


def crossproduct(x1, x2):
    a = x1[1]*x2[2] - x1[2]*x2[1]
    b = x1[2]*x2[0] - x1[0]*x2[2]
    c = x1[0]*x2[1] - x1[1]*x2[0]
    return [a, b, c]


def inverse(A):  # noqa
    det = (
        A[0][0]*A[1][1]*A[2][2] +
        A[1][0]*A[2][1]*A[0][2] +
        A[2][0]*A[0][1]*A[1][2] -
        A[0][2]*A[1][1]*A[2][0] -
        A[1][2]*A[2][1]*A[0][0] -
        A[2][2]*A[0][1]*A[1][0])
    if det:
        res = [
            [
                A[1][1]*A[2][2]/det - A[1][2]*A[2][1]/det,
                A[0][2]*A[2][1]/det - A[0][1]*A[2][2]/det,
                A[0][1]*A[1][2]/det - A[0][2]*A[1][1]/det
                ],
            [
                A[1][2]*A[2][0]/det - A[1][0]*A[2][2]/det,
                A[0][0]*A[2][2]/det - A[0][2]*A[2][0]/det,
                A[0][2]*A[1][0]/det - A[0][0]*A[1][2]/det
                ],
            [
                A[1][0]*A[2][1]/det - A[1][1]*A[2][0]/det,
                A[0][1]*A[2][0]/det - A[0][0]*A[2][1]/det,
                A[0][0]*A[1][1]/det - A[0][1]*A[1][0]/det
                ]
            ]
    else:
        res = None
    return res


def matplus(A, B):   # noqa
    if type(A) != list and type(B) != list:
        res = A+B
    elif type(A[0]) != list and type(B[0]) != list:
        res = [A[i] + B[i] for i in range(len(A))]
    else:
        try:
            res = [[A[i][j]+B[i][j] for j in range(len(B[0]))]
                   for i in range(len(A))]
        except TypeError:
            raise TypeError(
                "Matrices dimensions not compatible: %s and %s" % (A, B))
    return res


def matminus(A, B):  # noqa
    return matplus(A, matmul(-1, B))


def matmul(A, B):  # noqa
    if type(A) != list and type(B) != list:
        res = A*B
    elif type(A) != list:
        if type(B[0]) == list:
            res = [[A * B[i][j] for j in range(len(B[0]))]
                   for i in range(len(B))]
        else:
            res = [A * B[i] for i in range(len(B))]
    elif type(B) != list:
        res = matmul(B, A)
    elif type(A[0]) == list and type(B[0]) == list:
        res = [[sum(A[i][k]*B[k][j] for k in range(len(B)))
               for j in range(len(B[0]))] for i in range(len(A))]
    elif type(A[0]) == list:
        res = [sum(A[i][k]*B[k] for k in range(len(B)))
               for i in range(len(A))]
    elif type(B[0]) == list:
        res = [sum(A[k]*B[k][j] for k in range(len(B)))
               for j in range(len(B[0]))]
    else:
        res = dotproduct(A, B)
    return res


def matkronproduct(A, B):  # noqa
    return [
        [
            A[i][j]*B[k][l] for j in range(len(A[0])) for l in range(len(B[0]))
            ]
        for i in range(len(A)) for k in range(len(B))
        ]


def veckronproduct(x1, x2):
    return [[x1[i]*x2[j] for j in range(len(x2))] for i in range(len(x1))]


def transpose(A):  # noqa
    if type(A[0]) == list:
        res = [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]
    else:
        res = [[A[i]] for i in range(len(A))]
    return res


if __name__ == '__main__':
    A = [[1, 2, 3], [4, 5, 6], [7, 8, 10]]
    B = [[1, 2], [3, 4]]
    x1 = [1, 2, 3]
    x2 = [4, 5, 6]
    print A
    print transpose(A)
    print dotproduct(x1, x2)
    print crossproduct(x1, x2)
    print matkronproduct(A, B)
    print veckronproduct(x1, x2)
    print inverse(A)
    print matmul(A, x1)
    print matmul(x2, A)
    print transpose(x1)
    print matmul(3, x1)
    print matmul(A, 3)
