#!/usr/bin/python

# this script take .obj files in and add surface data to the .obj file

import argparse

def to_tetra(t_line):
    tetra = t_line.split(" ")
    tetra = list(lambda x: int(x), tetra[1:])
    assert( len(tetra) == 4 )
    return tetra

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="the .obj file to add surface to.")
    args = parser.parse_args()

    with fh = open(args.file, 'r'):
        s = fh.readlines()

        tetra = list( filter(lambda x: x.startswith("t"), s ) )
        face = list( filter(lambda x: x.startswith("f"), s ) )

        fh.close()

        if face:
            print("There is already faces in %s, please check before add surfaces." % args.file)
            return 0

        tetra = list( map(lambda x: )


