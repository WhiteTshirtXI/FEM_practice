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

def read_file(filename):

    with open(filename, 'r') as fh:
        s = fh.readlines();

        vertex = list( filter(lambda x: x.startswith("v"), s ) )
        face = list( filter(lambda x: x.startswith("f"), s ) )

        # transfer vertex line to vertex data
        print("transfer vertex...")
        vertex = list( map(to_vertex, vertex) )

        # transfer face line to face data
        print("transfer face...")
        face = list( map(to_face, face) )

        fh.close()

        return vertex, face

def write_hole(h, h_id, fh):
    for i in h:
        fh.write("h %d %d %d %d\n" % (h_id, i[0], i[1], i[2]))

vertices = []
faces = []
face_color = []
vf_list = []   # pair of (vertex, face)
v_end = []   # end location in vf_list of vertex faces


def fludfill(face_id, color):

    # waiting list of face_ids
    waiting_list = [face_id]

    face_color[face_id] = color;

    while waiting_list:
        f_expand = waiting_list.pop()

        # find nearby
        f_list = []
        for v in faces[f_expand]:
            for v_id, f_id in vf_list[ v_end[v-1] : v_end[v] ]:
                f_list.append(f_id)

        # color those uncolored
        f_list = list( filter( lambda f_id: face_color[f_id] == -1, f_list))
        for f_id in f_list:
            face_color[f_id] = color

        # add to quene
        waiting_list += f_list


vertex_boudary = []



def face_is_boundary(face):
    if not vertex_boudary[face[0] - 1]:
        return False

    if not vertex_boudary[face[1] - 1]:
        return False

    if not vertex_boudary[face[2] - 1]:
        return False

    return True


if __name__=="__main__":

    filename = "extrude_FDM.vtf"
    def vertex_is_boundary(v):
        return v[0] == 0

    vertices, faces = read_file(filename)

    print("building vertex-face list")
    for i, f in enumerate(faces):
        vf_list += [(f[0], i), (f[1], i), (f[2], i)]

    print("total %d v-f pairs, sorting ..." % len(vf_list))
    vf_list.sort(key = lambda x: x[0])

    print("building vertex_end list")
    count = 0
    v_end = []
    v_old_id = 0
    for v_id, f_id in vf_list:
        if v_old_id < v_id:
                        while (v_old_id < v_id):
                            v_old_id += 1
                            v_end.append(count)
        count += 1
    while (v_old_id < len(vertices)+1):
        v_old_id += 1
        v_end.append(count)
    assert(len(v_end) == len(vertices)+1), "v_end %d, vertices %d" % (len(v_end), len(vertices))

    face_color = [ -1 ] * len(faces)
    # color those unneeded as 0
    vertex_boudary = list(map(vertex_is_boundary, vertices))
    for i in range(len(faces)):
        if face_is_boundary(faces[i]):
            face_color[i] = 0

    # fludfill
    seed = 1
    uncolored = list(filter(lambda x: face_color[x] == -1, range(len(faces))))
    while uncolored:
        print("fludfill for color %d" % seed)
        face_id = uncolored.pop()
        fludfill(face_id, seed)
        seed += 1
        uncolored = list(filter(lambda x: face_color[x] == -1, range(len(faces))))

    # output
    face_number = max(face_color)
    for color in range(1, seed):
        out_faces = list( map(
                    lambda x: faces[x],
                    filter(lambda x: face_color[x] == color, range(len(faces)))
                    ))
        wh = open("hole_%d" % color, "w")
        write_hole(out_faces, color, wh)
        wh.close()
