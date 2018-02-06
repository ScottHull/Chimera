from Chimera import mesh, box

m = mesh.Mesh(x=3, y=3, z=3, spatial_res=0.2, verbose=True).build()
b = box.Box(mesh=m)