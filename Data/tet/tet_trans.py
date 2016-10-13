#!/usr/bin/python

import re
import argparse

ELEMENT_BIAS = 1

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help = "source *.tet file path")

    args = parser.parse_args()

    objname = re.split("\W", args.filename)[-2]
    print(objname)

    if (".tet" not in args.filename):
        print("file %s is not .tet file.\n" % args.filename)

    fh = open(args.filename, 'r')

    n_v = int( fh.readline().split(" ")[0] );
    t_v = int( fh.readline().split(" ")[0] );

    wh = open( objname + ".node", 'w')
    wh.write( "%d 3 0 0\n" % n_v )
    for i in range(n_v):
        pos = list(map(
            lambda x: float(x),
            fh.readline().split(" ")
            ))
        wh.write("%d %f %f %f\n" % ( i+1, *pos ) )
    wh.close()

    wh = open( objname + ".ele", 'w')
    wh.write( "%d 4 0\n" % t_v )
    for i in range(t_v):
        ele = list(map(
            lambda x: int(x) + ELEMENT_BIAS,
            fh.readline()[2:].split(" ")
            ))
        wh.write("%d %d %d %d %d\n" % (i+1, *ele) )
    wh.close()

    fh.close()

