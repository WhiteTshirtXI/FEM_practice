import re
import argparse
import numpy as np

# return v position
def to_vertex(v_line):
    vertex = v_line.split(" ")
    vertex = list( map(lambda x: float(x), vertex[1:]))
    assert( len(vertex) == 3 )
    return np.array(vertex)

# return surface pair (set(face_v), face_v)
def to_face(f_line):
    face = f_line.split(" ")
    face = list( map(lambda x: int(x), face[1:]))
    assert( len(face) == 3 )
    return face

# return tetra arange
def to_tetra(t_line):
    tetra = t_line.split(" ")
    tetra = list(map(lambda x: int(x), tetra[1:]))
    assert( len(tetra) == 4 )
    return tetra

def read_file(filename):

    with open(filename, 'r') as fh:
        s = fh.readlines();

        vertex = list( filter(lambda x: x.startswith("v"), s ) )
        face = list( filter(lambda x: x.startswith("f"), s ) )
        tetra = list( filter(lambda x: x.startswith("t"), s ) )

        # transfer vertex line to vertex data
        print("transfer vertex...")
        vertex = list( map(to_vertex, vertex) )

        # transfer face line to face data
        print("transfer face...")
        face = list( map(to_face, face) )

        # transfer tetra line to tetra data and arange the 4 vertices in order
        print("transfer tetra...")
        tetra = list( map(to_tetra, tetra) )

        fh.close()

        return vertex, face, tetra

# wirte fixed vertices to filehandle
def write_fixed(x, fh):
    fh.write("x")
    for i in x:
        fh.write(" %d" % i)
    fh.write("\n")

# wirte rigid bodies to filehandle
def write_rigid(r, fh):
    fh.write("r")
    for i in r:
        fh.write(" %d" % i)
    fh.write("\n")

# wirte hole to filehanle, each hole is list of faces
def write_hole(h, h_id, fh):
    for i in h:
        fh.write("h %d %d %d %d\n" % (h_id, i[0], i[1], i[2]))

if __name__=="__main__":

    filename = "penbox.vtf"

    # exame file
    with open(filename, 'r') as fh:
        s = fh.readlines();
        x = list( filter( lambda x: x.startswith("x"), s) )
        r = list( filter( lambda x: x.startswith("r"), s) )
        h = list( filter( lambda x: x.startswith("h"), s) )

        if x:
            print("There already exists fixed vertices in %s" % filename)
            return 0
        if r:
            print("There already exists rigid bodies in %s" % filename)
            return 0
        if h:
            print("There already exists hole information in %s" % filename)
            return 0

    # read in vertices, surface and tetra data
    v, f, t = read_file(filename)

    # some operations
    fixed_ids = []
    for i in range(len(v)):
        # y value equals -1
        if v[i][1] == -1:
            fixed_ids.append(i+1)

    rigid = []
    for i in range(len(v)):
        # y value equals 1
        if v[i][1] == 1:
            rigid.append(i+1)

    hole = []
    # detect inner wall surface vertices

    def inbox(pos, x_lim, y_lim, z_lim):
        if (pos[0] < x_lim[0]) or (pos[0] > x_lim[1]):
            return False
        if (pos[1] < y_lim[0]) or (pos[1] > y_lim[1]):
            return False
        if (pos[2] < z_lim[0]) or (pos[2] > z_lim[1]):
            return False
        return True

    inner = list( map( lambda x: inbox(x, [-0.9, 0.9], [-1.0, 0.9], [-0.9, 0.9]), v) )

    # detect hole face when given vertices inner
    def inset(x, v_set):
        return ( v_set[x[0] - 1] and  v_set[x[1] - 1] and  v_set[x[2] - 1] )

    hole = list( filter( lambda x: inset(x, inner), f) )

    with open(filename, 'a') as wh:
        write_fixed(fixed_ids, wh)
        write_rigid(rigid, wh)
        write_hole(hole, 1, wh)
