#!/usr/bin/env python
import math
import sys

__author__ = 'uemoto'

def dot(p, q):
    return p[0]*q[0]+p[1]*q[1]+p[2]+q[2]

def cross(p, q):
    return [p[1]*q[2]-p[2]*q[1],
            p[2]*q[0]-p[0]*q[2],
            p[0]*q[1]-p[1]*q[0]]

def det(A):
    return (+A[0][0]*A[1][1]*A[2][2]
            -A[0][0]*A[1][2]*A[2][1]
            +A[0][1]*A[1][2]*A[2][0]
            -A[0][1]*A[1][0]*A[2][2]
            +A[0][2]*A[1][0]*A[2][1]
            -A[0][2]*A[1][1]*A[2][0])

def inv(A):
    inv_d = 1.0 / det(A)
    i00 = +(A[1][1]*A[2][2]-A[1][2]*A[2][1])*inv_d
    i10 = -(A[1][0]*A[2][2]-A[1][2]*A[2][0])*inv_d
    i20 = +(A[1][0]*A[2][1]-A[1][1]*A[2][0])*inv_d
    i01 = -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*inv_d
    i11 = +(A[0][0]*A[2][2]-A[0][2]*A[2][0])*inv_d
    i21 = -(A[0][0]*A[2][1]-A[0][1]*A[2][0])*inv_d
    i02 = +(A[0][1]*A[1][2]-A[0][2]*A[1][1])*inv_d
    i12 = -(A[0][0]*A[1][2]-A[0][2]*A[1][0])*inv_d
    i22 = +(A[0][0]*A[1][1]-A[0][1]*A[1][0])*inv_d
    return [[i00, i01, i02], [i10, i11, i12], [i20, i21, i22]]

def prod(A, p):
    ap0 = dot(A[0], p)
    ap1 = dot(A[1], p)
    ap2 = dot(A[2], p)
    return [ap0, ap1, ap2]

class Chgcar:
    # Class for the data-storage of CHGCAR
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
        # Scalar parameters
        aa = dot(self.a, self.a)
        ab = dot(self.a, self.b)
        ac = dot(self.a, self.c)
        bb = dot(self.b, self.b)
        bc = dot(self.b, self.c)
        cc = dot(self.c, self.c)
        self.matrix = [[aa, ab, ac], [ab, bb, bc], [ac, bc, cc]]
        self.imatrix = inv(self.matrix)

    def get_data(self, s, ix, iy, iz):
        ix = ix % self.nx
        iy = iy % self.ny
        iz = iz % self.nz
        return self.rho[s][ix + self.nx * (iy + self.ny * iz)]


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
        for (jx, px) in enumerate([1.0-dx, dx]):
            for (jy, py) in enumerate([1.0-dy, dy]):
                for (jz, pz) in enumerate([1.0-dz, dz]):
                    w = px*py*pz
                    f = get_data(s, ix+jx, iy+jy, iz+jz)
                    r+= w*f
        return result


    def get_position(self, x, y, z):
        r = [x, y, z]
        ar = dot(self.a, r)
        br = dot(self.b, r)
        cr = dot(self.c, r)
        return prod(self.imatrix, [ar, br, cr])
    s


# Standard colorbar
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
            r = 1.0
            g = 0
            b = 0
    return (r, g, b)






def main():
    with open("CHGCAR") as fh:
        chgcar = Chgcar(fh)
        sq1 = 0.0
        sq2 = 0.0
        num = 0
        with open("output.pam", "wb") as img:
            img.write("P7\n")
            img.write("WIDTH 300\n")
            img.write("HEIGHT 300\n")
            img.write("DEPTH 4\n")
            img.write("MAXVAL 255\n")
            img.write("TUPLTYPE RGB_ALPHA\n")
            img.write("ENDHDR\n")
            for i in xrange(300):
                for j in xrange(300):
                    x = i*0.01
                    y = j*0.01
                    z = 0.5
                    q = (chgcar.density(0, x, y, z))
                    sq1 += q
                    sq2 += q*q
                    num += 1
                    (r, g, b) = colorbar(q/10)
                    img.write(chr(int(r*255)))
                    img.write(chr(int(g*255)))
                    img.write(chr(int(b*255)))
                    img.write(chr(255))
            print "# Total %d, Average %e, Sigma %e" % (
                num, sq1 / num, math.sqrt((sq1*sq1-sq2))/num
            )
if __name__ == "__main__":
    main()
