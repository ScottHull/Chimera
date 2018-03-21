import numpy as np
import time
from Chimera.Chimera_3D import console


class Mesh:

    def __init__(self, spatial_res, x, y, z=None, verbose=True):
        self.verbose = verbose
        self.spatial_res = spatial_res
        if '.' in str(self.spatial_res):
            self.spatial_sigfigs = len(str(self.spatial_res)) - 1
        else:
            self.spatial_sigfigs = len(str(self.spatial_res))
        self.x = x
        self.y = y
        self.x_coords = None
        self.y_coords = None
        if z is not None:
            self.z = z
            self.z_coords = None
            self.dimensions = 3
        else:
            self.dimensions = 2



    def build(self, df=None):
        console.event("Generating model nodes...", verbose=self.verbose)
        t_start = time.time()
        x_coords = np.arange(0, self.x + self.spatial_res, self.spatial_res)
        y_coords = np.arange(0, self.y + self.spatial_res, self.spatial_res)
        z_vals = []
        nodes = []
        if self.dimensions is 3:
            console.nominal("Detected a 3-D system! Generating a {} x {} x {} model!".format(
                round(self.x / self.spatial_res, self.spatial_sigfigs),
                round(self.y / self.spatial_res, self.spatial_sigfigs),
                round(self.z / self.spatial_res, self.spatial_sigfigs)),
                verbose=self.verbose)
            z_coords = np.arange(0, self.z + self.spatial_res, self.spatial_res)
            for x in x_coords:
                for y in y_coords:
                    for z in z_coords:
                        x = round(x, self.spatial_sigfigs)
                        y = round(y, self.spatial_sigfigs)
                        z = round(z, self.spatial_sigfigs)
                        z_vals.append(z)
                        node = (x, y, z)
                        nodes.append(node)
        else:
            console.nominal("Detected a 2-D system! Generating a {} x {} model!".format(
                round(self.x / self.spatial_res, self.spatial_sigfigs),
                round(self.y / self.spatial_res, self.spatial_sigfigs)),
                verbose=self.verbose)
            for x in x_coords:
                for y in y_coords:
                    x = round(x, self.spatial_sigfigs)
                    y = round(y, self.spatial_sigfigs)
                    z_vals.append(y)
                    node = (x, y)
                    nodes.append(node)
        if df is not None:
            df['coords'] = nodes
            df['z_vals'] = z_vals
        time_task = time.time() - t_start
        console.event("Finished generating model vertices! (task took {}s)".format(time_task), verbose=self.verbose)
        return nodes


    def get_spatial_sigfigs(self):
        return self.spatial_sigfigs
