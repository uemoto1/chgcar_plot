#!/usr/bin/env python
import math
import sys

__author__ = 'uemoto'


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
        ab = [
            self.a[1]*self.b[2]-self.a[2]*self.b[1],
            self.a[2]*self.b[0]-self.a[0]*self.b[2],
            self.a[0]*self.b[1]-self.a[1]*self.b[0]
        ]
        abc = ab[0]*self.c[0]+ab[1]*self.c[1]+ab[2]*self.c[2]
        self.volume = abc


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
        Dx = 1.0 - dx
        Dy = 1.0 - dy
        Dz = 1.0 - dz
        r000 = self.get_data(s, ix, iy, iz)
        r001 = self.get_data(s, ix, iy, iz+1)
        r010 = self.get_data(s, ix, iy+1, iz)
        r011 = self.get_data(s, ix, iy+1, iz+1)
        r100 = self.get_data(s, ix+1, iy, iz)
        r101 = self.get_data(s, ix+1, iy, iz+1)
        r110 = self.get_data(s, ix+1, iy+1, iz)
        r111 = self.get_data(s, ix+1, iy+1, iz+1)
        return (r000*Dx*Dy*Dz + r001*Dx*Dy*dz + r010*Dx*dy*Dz + r011*Dx*dy*dz
            + r100*dx*Dy*Dz + r101*dx*Dy*dz + r110*dx*dy*Dz + r111*dx*dy*dz)




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
