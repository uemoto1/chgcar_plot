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
        if (0 <= dx <= 1) and (0 <= dy <= 1) and (0 <= dz <= 1):
            pass
        else:
            print dx, dy, dz
            sys.exit(0)
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


def main():
    with open("CHGCAR") as fh:
        chgcar = Chgcar(fh)
        for i in xrange(300):
            for j in xrange(300):
                x = i*0.01
                y = j*0.01
                z = 0.5
                print "%f %f %e" % (x, y, chgcar.density(0, x, y, z))
        for (x, y, z) in chgcar.atom:
            print "# %f %f %f" % (x, y, z)

if __name__ == "__main__":
    main()
