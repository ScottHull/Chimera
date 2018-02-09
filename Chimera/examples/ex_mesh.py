from Chimera import box

b = box.Box(evolution_time=20, conduction=True, settling_mode='stokes terminal', radioactivity=True,
            chemistry=True, verbose=True)
c = b.build(x=5, y=5, z=5, spatial_res=0.1)
