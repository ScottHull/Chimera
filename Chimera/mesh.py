import numpy as np
import time
from . import console


# This is currently for use in a 1D model...



class Mesh:

    def __init__(self, spatial_res, z, verbose=True):
        self.verbose = verbose
        self.spatial_res = spatial_res
        if '.' in str(self.spatial_res):
            self.spatial_sigfigs = len(str(self.spatial_res)) - 1
        else:
            self.spatial_sigfigs = len(str(self.spatial_res))
        self.z = z
        self.z_coords = None
        self.dimensions = 1


    def build_linear(self, df=None):
        console.event("Generating model nodes...", verbose=self.verbose)
        t_start = time.time()
        nodes = []
        z_coords = np.arange(0.0, self.z + self.spatial_res, self.spatial_res)
        console.event("Detected a 1D model.  Generating {} nodes...".format(len(z_coords)), verbose=self.verbose)
        for k in z_coords:
            nodes.append(round(k, self.spatial_sigfigs))
        if df is not None:
            df['coords'] = nodes
        console.event("Finished generating model nodes! (task took {}s)".format(
            t_start - time.time()), verbose=self.verbose)
        return nodes

    def get_spatial_sigfigs(self):
        return self.spatial_sigfigs


