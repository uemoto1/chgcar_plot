#!/usr/bin/env python
from optparse import OptionParser
from numpy import array




__author__ = 'uemoto'





def read_chgcar(chgcar):
    flag_keyword = False
    flag_data = False
    flag_skip = False
    data = []
    temp = []
    for (i, buff) in enumerate(chgcar):
        line = buff.strip()
        if flag_skip:
            print line
            if keyword in line:
                flag_skip = False
        elif flag_data:
            token = line.split()
            temp += map(float, token)
            print len(temp), nx, ny, nz, nsize
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
            flag_data = True
        elif len(line) == 0:
            flag_keyword = True
    return data



def main()
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", type="string", default=None, help="CHGCAR file")
    opts, args = parser.parse_args()
    print opts, args




if __name__ == "__main__":
    main()