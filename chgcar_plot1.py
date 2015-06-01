#!/usr/bin/env python
import optparse
import math
import sys
import os

__author__ = 'uemoto'

# Dot product of 3-dimensional vectors


def dot(p, q):
    return p[0] * q[0] + p[1] * q[1] + p[2] * q[2]

# Cross product of 3-dimensional vectors


def cross(p, q):
    return [p[1] * q[2] - p[2] * q[1],
            p[2] * q[0] - p[0] * q[2],
            p[0] * q[1] - p[1] * q[0]]

# Normalize


def normal(p):
    ip = 1.0 / math.sqrt(dot(p, p))
    return [p[0] * ip, p[1] * ip, p[2] * ip]

# Determinant of 3x3 matrix


def det(A):
    return (+A[0][0] * A[1][1] * A[2][2]
            - A[0][0] * A[1][2] * A[2][1]
            + A[0][1] * A[1][2] * A[2][0]
            - A[0][1] * A[1][0] * A[2][2]
            + A[0][2] * A[1][0] * A[2][1]
            - A[0][2] * A[1][1] * A[2][0])

# Inverse matric of A


def inv(A):
    inv_det = 1.0 / det(A)
    i00 = +(A[1][1] * A[2][2] - A[1][2] * A[2][1]) * inv_det
    i01 = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) * inv_det
    i02 = +(A[0][1] * A[1][2] - A[0][2] * A[1][1]) * inv_det
    i10 = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) * inv_det
    i11 = +(A[0][0] * A[2][2] - A[0][2] * A[2][0]) * inv_det
    i12 = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) * inv_det
    i20 = +(A[1][0] * A[2][1] - A[1][1] * A[2][0]) * inv_det
    i21 = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) * inv_det
    i22 = +(A[0][0] * A[1][1] - A[0][1] * A[1][0]) * inv_det
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
        flag = 0  # 0:header, 1:position, 2:keyword, 3:charge, 4:skip
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
                flag = 1  # Reading Atom positon
            elif 8 <= i:
                if flag == 1:
                    if 0 < len(line):
                        pos = map(float, line.split())
                        self.atom.append(pos)
                    else:
                        temp = []  # Initialize buffer
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
        for (jx, px) in enumerate([1.0 - dx, dx]):
            for (jy, py) in enumerate([1.0 - dy, dy]):
                for (jz, pz) in enumerate([1.0 - dz, dz]):
                    w = px * py * pz
                    f = self.get_data(s, ix + jx, iy + jy, iz + jz)
                    r += w * f
        return r


def colorbar_mono(c):
    if c < 0.0:
        r = 0.25
        g = 0.25
        b = 0.25
    elif c < 1.00:
        r = 0.25 + c * 0.50
        g = 0.25 + c * 0.50
        b = 0.25 + c * 0.50
    else:
        r = 0.75
        g = 0.75
        b = 0.75
    return (r, g, b)


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


def plot(chg, r0, r1, r2, xr, yr, output="chgcar.pam", factor=0.0005, sample=0.1, mode=0, flag=False):
    [x0, y0, z0] = chg.abc2xyz(r0)
    d1 = chg.abc2xyz(r1)
    d2 = chg.abc2xyz(r2)
    u = normal(d1)
    w = normal(cross(d1, d2))
    v = normal(cross(w, u))
    x1, x2 = xr
    y1, y2 = yr
    nu = int((x2 - x1) / sample)
    nv = int((y2 - y1) / sample)
    print u, v, w
    with open(output, "wb") as out:
        out.write("P7\n")
        out.write("WIDTH %d\n" % nu)
        out.write("HEIGHT %d\n" % nv)
        out.write("DEPTH 4\n")
        out.write("MAXVAL 255\n")
        out.write("TUPLTYPE RGB_ALPHA\n")
        out.write("ENDHDR\n")
        for j in xrange(nv):
            vv = sample * j + y1
            for i in xrange(nu):
                uu = sample * i + x1
                x = x0 + u[0] * uu + v[0] * vv
                y = y0 + u[1] * uu + v[1] * vv
                z = z0 + u[2] * uu + v[2] * vv
                [a, b, c] = chg.xyz2abc([x, y, z])
                if (0 <= a < 1) and (0 <= b < 1) and (0 <= c < 1) or flag:
                    if mode == 0:
                        q = factor * chg.density(0, a, b, c)
                    else:
                        q = factor * (chg.density(mode, a, b, c)) + 0.5
                    # Plot the calculated results
                    [r, g, b] = colorbar(q)
                    out.write(chr(int(255 * r)))
                    out.write(chr(int(255 * g)))
                    out.write(chr(int(255 * b)))
                    out.write(chr(255))
                else:
                    out.write(chr(0))
                    out.write(chr(0))
                    out.write(chr(0))
                    out.write(chr(0))


def main():
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", dest="input", type=str,
                      default="CHGCAR", help="input 'CHGCAR' file")
    parser.add_option("-o", "--output", dest="output", type=str,
                      default="chgcar.pam", help="output 'pam' image file")
    parser.add_option("-c", "--center", dest="center", type=str,
                      default="0,0,0", help="center position")
    parser.add_option("-u", "--uvec", dest="uvec", type=str,
                      default="0,1,0", help="vector of u-direction in figure")
    parser.add_option("-v", "--vvec", dest="vvec", type=str,
                      default="0,0,1", help="vector of v-direction in figure")
    parser.add_option("-U", "--urange", dest="urange", type=str,
                      default="-1:10", help="range of u in plot")
    parser.add_option("-V", "--vrange", dest="vrange", type=str,
                      default="-1:16", help="range of v in plot")
    parser.add_option("-s", "--sample", dest="sample", type=float,
                      default=0.05, help="range of v in plot")
    parser.add_option("-m", "--mode", dest="mode", type=int,
                      default=0, help="0:Charge density, 1,2,3: Spin density of x,y,z")
    parser.add_option("-a", "--allarea", dest="allarea", action="store_true",
                      default=False, help="plot all region")
    (opts, args) = parser.parse_args()

    if not os.path.isfile(opts.input):
        sys.stderr.write("Error! '%s' is not found!\n" % (opts.input))
        sys.exit(-1)

    c0 = map(float, opts.center.split(","))
    uv = map(float, opts.uvec.split(","))
    vv = map(float, opts.vvec.split(","))
    ur = map(float, opts.urange.split(":"))
    vr = map(float, opts.vrange.split(":"))

    with open(opts.input) as fh:
        chgcar = Chgcar(fh)
        plot(
            chgcar, c0, uv, vv, ur, vr, output=opts.output, factor=0.001, sample=opts.sample,
            mode=opts.mode, flag=opts.allarea)


if __name__ == "__main__":
    main()
