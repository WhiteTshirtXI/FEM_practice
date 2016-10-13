import numpy as np
import matplotlib.pyplot as plt

def sepLine(line):
    l = line.split(" ")
    return list(map(lambda x: float(x), l))

if __name__=='__main__':
    data = open('cylinder.txt').readlines()

    d_list = []
    f_list = []

    for line in data:
        d, f = sepLine(line)
        d_list.append(d)
        f_list.append(f)

    d_list = np.array(d_list)
    f_list = np.array(f_list)

    A = np.vstack([d_list, np.ones(len(d_list))]).T
    m, c = np.linalg.lstsq(A, f_list)[0]

    print("Young's modulus: %f", m*4.0/np.pi)

    plt.plot(d_list, f_list, 'o')
    plt.xlabel("displacement / m")
    plt.ylabel("Force / N")
    plt.show()
