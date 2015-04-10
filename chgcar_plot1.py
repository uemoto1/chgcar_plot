#!/usr/bin/env python
from optparse import OptionParser
from numpy import array




__author__ = 'uemoto'





def read_chgcar(chgcar):
    flag_keyword = False
    flag_data = False
    flag_skip = False
    data = []
    for i, buff in enumerate(chgcar):
        line = buff.strip()
        if flag_skip:
            print line
            if keyword in line:
                flag_skip = False
        elif flag_data:
            token = line.split()
            temp += map(float, token)
            if nsize <= len(temp):
                flag_skip = True
                data.append(temp)
                temp = []
        elif flag_keyword:
            keyword = line
            token = keyword.split()
            nx = int(token[0])
            ny = int(token[1])
            nz = int(token[2])
            nsize = nx*ny*nz
            temp = []
            flag_data = True
        elif len(line) == 0:
            flag_keyword = True
        elif i == 2:
            vec_a = map(float, line.split())
        elif i == 3:
            vec_b = map(float, line.split())
        elif i == 4:
            vec_c = map(float, line.split())
    lattice = [vec_a, vec_b, vec_c]
    size = [nx, ny, nz]
    return [lattice, size, data]


def main()
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", type="string", default="CHGCAR", help="CHGCAR file")
    opts, args = parser.parse_args()
    print opts, args




if __name__ == "__main__":
    main()