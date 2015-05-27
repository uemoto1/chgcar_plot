#!/usr/bin/env python
__author__ = 'uemoto'

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
        ab = [
            self.a[1]*self.b[2]-self.a[2]*self.b[1],
            self.a[2]*self.b[0]-self.a[0]*self.b[2],
            self.a[0]*self.b[1]-self.a[1]*self.b[0]
        ]
        abc = ab[0]*self.c[0]+ab[1]*self.c[1]+ab[2]*self.c[2]
        self.volume = abc


    def get_value(s, i, j, k):
        


def main():
    with open("CHGCAR") as fh:
        chgcar = Chgcar(fh)



if __name__ == "__main__":
    main()
