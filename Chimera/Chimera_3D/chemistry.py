import pandas as pd
import numpy as np
from . import console, backends


class Chemistry:

    def __init__(self, box):
        self.box = box
        self.mesh = box.mesh
        self.objects = box.objects
        self.coords = np.array(self.mesh['coords'])
        self.verbose = box.verbose

    def insert_matrix_composition(self, material):
        pass

