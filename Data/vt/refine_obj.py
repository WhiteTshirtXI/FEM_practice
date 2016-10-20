#!/usr/bin/python

# this script take .obj files in and add rearange the vertex ids in tetra
# so that for tetra [a, b, c, d], the face [a, b, c] is towards outside.
# then generate surface data to the .obj file

import argparse
import re
import numpy as np

# return v position
def to_vertex(v_line):
    vertex = v_line.split(" ")
    vertex = list( map(lambda x: float(x), vertex[1:]))
    assert( len(vertex) == 3 )
    return np.array(vertex)

# return surface pair (set(face_v), face_v)
def to_surface(t):
    face = []
    face.append([t[0], t[1], t[2]])
    face.append([t[0], t[2], t[3]])
    face.append([t[0], t[3], t[1]])
    face.append([t[3], t[2], t[1]])

    surf = list( map( lambda x: (sorted(x), x), face) )
    return surf

# return tetra arange
def to_tetra(t_line):
    tetra = t_line.split(" ")
    tetra = list(map(lambda x: int(x), tetra[1:]))
    assert( len(tetra) == 4 )
    return tetra

# exame if the order of vertices in tetra is correct, then
# return the correct ordered tetra
def exame_and_switch(v_pos, tetra):
    # notice that the location in list = v_id - 1
    r = list( map(lambda x: v_pos[x - 1], tetra))

    if np.dot(np.cross(r[0] - r[3], r[1] -r[3]), r[2] - r[3]) < 0:
        return [tetra[1], tetra[0], tetra[2], tetra[3]]

    return tetra


# this function take .obj file in and rearange the vertex ids in tetra
# so that for tetra [a, b, c, d], the face [a, b, c] is towards outside
def refine_tetra(in_file, out_file):
    with open(in_file, 'r') as fh:
        s = fh.readlines()

        vertex = list( filter(lambda x: x.startswith("v"), s ) )
        tetra = list( filter(lambda x: x.startswith("t"), s ) )

        # transfer vertex line to vertex data
        print("transfer vertex...")
        v_pos = list( map(to_vertex, vertex) )

        # transfer tetra line to tetra data and arange the 4 vertices in order
        print("transfer tetra...")
        t_aranged = list( map(lambda x: exame_and_switch(v_pos, to_tetra(x)), tetra) )

        print("writing...")
        wh = open(out_file, "w")

        for v_line in vertex:
            wh.write(v_line)

        for tetra in t_aranged:
            wh.write("t %d %d %d %d\n" % (tetra[0], tetra[1], tetra[2], tetra[3]))

        wh.close()

# generate surface for aranged tetra file .obj
def add_surface(in_file):

    with open(in_file, 'r') as fh:
        s = fh.readlines()

        tetra = list( filter(lambda x: x.startswith("t"), s ) )
        face = list( filter(lambda x: x.startswith("f"), s ) )

        fh.close()

        if face:
            print("There is already faces in %s, please check before add surfaces." % in_file)
            return 0

        tetra = list( map(to_tetra, tetra) )

        print('building tetra faces...')
        surf = []
        for t in tetra:
            surf += to_surface(t)
        print("toal %d faces created" % len(surf))

        print("sorting face pairs...")
        surf.sort(key=lambda x: x[0])

        print("picking single faces...")
        single = []
        while len(surf)>1:
            if (surf[-1][0] == surf[-2][0]):
                surf.pop()
                surf.pop()
            else:
                single.append(surf.pop()[1])
        if surf:
            single.append(surf.pop()[1])

        print("total %d faces of surface" % len(single))

        print("writing...")
        wh = open(in_file, 'a')
        for f in single:
            wh.write("f %d %d %d\n" % (f[0], f[1], f[2]))
        wh.close()



if __name__=="__main__":

    # test exam_and_switch
    a = np.array([0,0,0])
    b = np.array([1,0,0])
    c = np.array([0,1,0])
    d = np.array([0,0,1])
    v_pos = [a,b,c,d]
    tetra = [1, 2, 3, 4]
    assert( exame_and_switch(v_pos, tetra) == [2, 1, 3,4])


    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="the .vt file to add surface to.")
    args = parser.parse_args()

    filename = re.split("\W", args.file)[-2]
    print("writing to %s.obj" % filename)
    refine_tetra(args.file, filename + '.obj')
    add_surface(filename + '.obj')
