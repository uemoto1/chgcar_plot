#!/usr/bin/env python
import math
import sys

__author__ = 'uemoto'

# Dot product of 3-dimensional vectors
def dot(p, q):
    return p[0]*q[0]+p[1]*q[1]+p[2]*q[2]

# Cross product of 3-dimensional vectors
def cross(p, q):
    return [p[1]*q[2]-p[2]*q[1],
            p[2]*q[0]-p[0]*q[2],
            p[0]*q[1]-p[1]*q[0]]

# Normalize
def normal(p):
    ip = 1.0 / math.sqrt(dot(p, p))
    return [p[0]*ip, p[1]*ip, p[2]*ip]

# Determinant of 3x3 matrix
def det(A):
    return (+A[0][0]*A[1][1]*A[2][2]
            -A[0][0]*A[1][2]*A[2][1]
            +A[0][1]*A[1][2]*A[2][0]
            -A[0][1]*A[1][0]*A[2][2]
            +A[0][2]*A[1][0]*A[2][1]
            -A[0][2]*A[1][1]*A[2][0])

# Inverse matric of A
def inv(A):
    inv_det = 1.0 / det(A)
    i00 = +(A[1][1]*A[2][2]-A[1][2]*A[2][1])*inv_det
    i01 = -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*inv_det
    i02 = +(A[0][1]*A[1][2]-A[0][2]*A[1][1])*inv_det
    i10 = -(A[1][0]*A[2][2]-A[1][2]*A[2][0])*inv_det
    i11 = +(A[0][0]*A[2][2]-A[0][2]*A[2][0])*inv_det
    i12 = -(A[0][0]*A[1][2]-A[0][2]*A[1][0])*inv_det
    i20 = +(A[1][0]*A[2][1]-A[1][1]*A[2][0])*inv_det
    i21 = -(A[0][0]*A[2][1]-A[0][1]*A[2][0])*inv_det
    i22 = +(A[0][0]*A[1][1]-A[0][1]*A[1][0])*inv_det
    return [[i00, i01, i02], [i10, i11, i12], [i20, i21, i22]]

# Matrix-vector production
def prod(A, p):
    ap0 = dot(A[0], p)
    ap1 = dot(A[1], p)
    ap2 = dot(A[2], p)
    return [ap0, ap1, ap2]


# Class for the data-storage of CHGCAR
class Chgcar:
    def __init__(self, chgcar_file):
        skip_keyword = "augmentation occupancies"
        # Read numerical data from CHGCAR format
        self.rho = []
        self.atom = []
        flag = 0 # 0:header, 1:position, 2:keyword, 3:charge, 4:skip
        for i, buff in enumerate(chgcar_file):
            line = buff.strip()
            if i == 0:
                self.title = line
            elif i == 1:
                self.scale = float(line)
            elif i == 2:
                self.a = map(float, line.split())
            elif i == 3:
                self.b = map(float, line.split())
            elif i == 4:
                self.c = map(float, line.split())
            elif i == 5:
                self.elem = line.split()
            elif i == 6:
                self.nelem = map(int, line.split())
            elif i == 7:
                self.mode = line
                flag = 1 # Reading Atom positon
            elif 8 <= i:
                if flag == 1:
                    if 0 < len(line):
                        pos = map(float, line.split())
                        self.atom.append(pos)
                    else:
                        temp = [] # Initialize buffer
                        flag = 2
                elif flag == 2:
                    start_keyword = line
                    size = map(int, line.split())
                    self.nx = size[0]
                    self.ny = size[1]
                    self.nz = size[2]
                    flag = 3
                elif flag == 3:
                    if skip_keyword in line:
                        self.rho.append(temp)
                        temp = []
                        flag = 4
                    else:
                        temp += map(float, line.split())
                elif flag == 4:
                    if start_keyword in line:
                        flag = 3
        # Lattice parameters
        self.volume = dot(self.a, cross(self.b, self.c))
        # Transform matrix
        self.matrix = [
            [self.a[0], self.b[0], self.c[0]],
            [self.a[1], self.b[1], self.c[1]],
            [self.a[2], self.b[2], self.c[2]],
        ]
        self.imatrix = inv(self.matrix)
        print self.imatrix

    def get_data(self, s, ix, iy, iz):
        ix = ix % self.nx
        iy = iy % self.ny
        iz = iz % self.nz
        return self.rho[s][ix + self.nx * (iy + self.ny * iz)]

    def xyz2abc(self, xyz):
        return prod(self.imatrix, xyz)

    def abc2xyz(self, abc):
        return prod(self.matrix, abc)

    def density(self, s, a, b, c):
        x = a * self.nx
        y = b * self.ny
        z = c * self.nz
        ix = int(math.floor(x))
        iy = int(math.floor(y))
        iz = int(math.floor(z))
        dx = x - ix
        dy = y - iy
        dz = z - iz
        r = 0.0
        # Cubic-linear interpolation
        for (jx, px) in enumerate([1.0-dx, dx]):
            for (jy, py) in enumerate([1.0-dy, dy]):
                for (jz, pz) in enumerate([1.0-dz, dz]):
                    w = px*py*pz
                    f = self.get_data(s, ix+jx, iy+jy, iz+jz)
                    r+= w*f
        return r





def colorbar(c):
    if c < 0.0:
        r = 0
        g = 0
        b = 0.5
    elif c < 0.125:
        r = 0
        g = 0
        b = 0.5 + (c - 0.0) * 4.0
    elif c < 0.375:
        r = 0.0
        g = (c - 0.125) * 4.0
        b = 1.0
    elif c < 0.625:
        r = (c - 0.375) * 4.0
        g = 1.0
        b = 1.0 - (c - 0.375) * 4.0
    elif c < 0.875:
        r = 1.0
        g = 1.0 - (c - 0.625) * 4.0
        b = 0.0
    elif c < 1.0:
        r = 1.0 - (c - 0.875) * 4.0
        g = 0
        b = 0
    else:
        r = 0.5
        g = 0
        b = 0
    return (r, g, b)


def plot(chg, r0, r1, r2, xr, yr, s=0.1):
    [x0, y0, z0] = chg.abc2xyz(r0)
    d1 = chg.abc2xyz(r1)
    d2 = chg.abc2xyz(r2)
    u = normal(d1)
    w = normal(cross(d1, d2))
    v = normal(cross(w, u))
    x1, x2 = xr
    y1, y2 = yr
    nu = int((x2-x1)/s)
    nv = int((y2-y1)/s)
    print u, v, w
    with open("output.pam", "wb") as out:
        out.write("P7\n")
        out.write("WIDTH %d\n" % nu)
        out.write("HEIGHT %d\n" % nv)
        out.write("DEPTH 4\n")
        out.write("MAXVAL 255\n")
        out.write("TUPLTYPE RGB_ALPHA\n")
        out.write("ENDHDR\n")
        for j in xrange(nv):
            for i in xrange(nu):
                x = x0+u[0]*s*i+v[0]*s*i
                y = y0+u[1]*s*i+v[1]*s*j
                z = z0+u[2]*s*i+v[2]*s*j
                [a, b, c] = chg.xyz2abc([x, y, z])
                if (0 <= a <= 1) and (0 <= b <= 1) and (0 <= c <= 1):
                    q = chg.density(0, a, b, c)
                    [r, g, b] = colorbar(q/1000)
                    out.write(chr(int(255*r)))
                    out.write(chr(int(255*g)))
                    out.write(chr(int(255*b)))
                    out.write(chr(255))
                else:
                    out.write(chr(0))
                    out.write(chr(0))
                    out.write(chr(0))
                    out.write(chr(0))




def main():
    with open("CHGCAR") as fh:
        chgcar = Chgcar(fh)
        plot(chgcar, [0.5, -0.5, -0.5],
            [0., 1., 0.], [0., 0., 1.],
            20.00, 20.00, 0.02)

if __name__ == "__main__":
    main()
