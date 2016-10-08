#!/usr/bin/python
# author: Like Ma, <milkpku>at<gmail>dot<com>

# transfer data from stellar to .obj file
# v for vertex
# t for tetrahedron

import argparse

def read_vertex(fh):
    vertex = []

    # get size
    v_size = fh.readline()
    v_size = list(filter(lambda x: x, v_size.split(' ')))
    v_size = int(v_size[0])

    for l in fh:

        print(l)

        if l.startswith("#"):
            continue

        v = l.split(' ')
        v = list(filter(lambda x: x, v))
        print(v)
        v = v[1:4]
        vertex.append(v)

    return v_size, vertex

def read_tetra(fh):
    tetra = []

    # get size
    t_size = fh.readline()
    t_size = list(filter(lambda x: x, t_size.split(' ')))
    t_size = int(t_size[0])

    for l in fh:
        if l.startswith("#"):
            continue

        t = l.split(' ')
        t = list(filter(lambda x: x, t))
        t = t[1:5]
        tetra.append(t)

    return t_size, tetra



if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="the prefix of .node file and .ele file")
    args = parser.parse_args()

    v_size, vertex = read_vertex( open(args.file + '.node', 'r') )
    t_size, tetra = read_tetra( open(args.file + '.ele', 'r') )

    print("vertex number: %d" % len(vertex))
    print("tetra number: %d" % len(tetra))

    fh = open(args.file + '.obj', 'w')

    for v in vertex:
        fh.write('v %s %s %s\n' % (v[0], v[1], v[2]))

    for t in tetra:
        fh.write('t %s %s %s %s\n' % (t[0], t[1], t[2], t[3]))

    fh.close()

