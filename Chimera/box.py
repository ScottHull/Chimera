import os
import pandas as pd
from . import console

class Box:

    def __init__(self, mesh, evolution_time, conduction=True, settling_mode='stokes terminal',
                 radioactivity='active', chemistry='active'):
        self.space = pd.DataFrame({
            'verticies': (i for i in mesh)
        })
        self.conduction = conduction
        self.settling_mode = settling_mode
        self.radioactivity = radioactivity
        self.chemistry = chemistry
        self.initial_time = evolution_time
        self.evolution_time = evolution_time

    def update(self):
        pass
