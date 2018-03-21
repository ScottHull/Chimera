import numpy as np

x = 5
y = 5
z = 5
spatial_res = 0.2
spatial_sigfigs = 2

def predict(node, nodes):
    i = node[0]
    j = node[1]
    k = node[2]
    max_k = (z / spatial_res) + 1
    max_j = ((y / spatial_res) + 1) * max_k

    predicted = int(((i / spatial_res) * max_j) + ((j / spatial_res) * max_k) + k/spatial_res)
    predicted_2 = round((((i*y*z)/(spatial_res**3)) + (((i*y) + (i*z) + (j*z))/(spatial_res**2)) + ((i+j)/spatial_res) + (k/spatial_res)))


    print(predicted, predicted_2)
    actual = nodes.index(node)
    print("{}x{}x{} = {} box".format(x, y, z, len(nodes)))
    print(predicted, actual)
    print(nodes[predicted], nodes[actual])

def build():
    x_coords = np.arange(0, x + spatial_res, spatial_res)
    y_coords = np.arange(0, y + spatial_res, spatial_res)
    z_coords = np.arange(0, z + spatial_res, spatial_res)
    nodes = []
    for i in x_coords:
        for j in y_coords:
            for k in z_coords:
                p = round(i, spatial_sigfigs)
                r = round(j, spatial_sigfigs)
                q = round(k, spatial_sigfigs)
                node = (p, r, q)
                nodes.append(node)
    return nodes


predict(node=(1.0, 0.0, 0.2), nodes=build())

