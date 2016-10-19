#/usr/bin/python

# this script take .obj file in and rearange the vertex ids in tetra
# so that for tetra [a, b, c, d], the face [a, b, c] is towards outside

import argparse
import numpy as np

# return tetra arange
def to_tetra(t_line):
    tetra = t_line.split(" ")
    tetra = list(map(lambda x: int(x), tetra[1:]))
    assert( len(tetra) == 4 )
    return tetra

# return v position
def to_vertex(v_line):
    vertex = v_line.split(" ")
    vertex = list( map(lambda x: float(x), tetra[1:]))
    assert( len(vertex) == 3 )
    return np.array(vertex)

def exame_and_switch(v_pos, tetra):
    r = list( map(lambda x: v_pos[x - 1], tetra))

    if np.dot(np.cross(r[0] - r[3], r[1] -r[3]), r[2] - r[3]) < 0:
        return [tetra[1], tetra[0], tetra[2], tetra[3]]

    return tetra

def refine_tetra(in_file, out_file):
    with open(in_file, 'r') as fh:
        s = fh.readlines()

        vertex = list( filter(lambda x: x.startswith("v"), s ) )
        tetra = list( filter(lambda x: x.startswith("t"), s ) )

        print("transfer vertex")
        v_pos = list( map(to_vertex, vertex) )

        print("transfer tetra")
        t_aranged = list( map(lambda x: exame_and_switch(v_pos, to_tetra(x)), tetra) )

        print("writing")
        wh = open(out_file, "w")

        for v_line in vertex:
            wh.write(v_line)

        for tetra in t_aranged:
            wh.write("t %d %d %d %d\n" % (tetra[0], tetra[1], tetra[2], tetra[3]))

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

    # begin transfer
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="the .obj file to add surface to.")
    args = parser.parse_args()

    refine_tetra(args.file, "output.obj")




