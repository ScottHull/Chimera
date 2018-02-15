from Chimera import box

x = 50
y = 50
z = 50
spatial_res = 0.5

b = box.Box(evolution_time=20, conduction=True, settling_mode='stokes terminal', radioactivity=True,
            chemistry=True, verbose=True)
b.build(spatial_res=spatial_res, x=x, y=y, z=z)
b.insert_matrix(material='test_matrix', temperature=2000, depth_range=[0.5, 49.5])
b.insert_object(material='test_object', temperature=2000, radius=0.01, x=1, y=1, z=1)
b.insert_boundary(temperature=3000, depth_range=[0.0, 0.5])
b.insert_boundary(temperature=4000, depth_range=[49.5, 50.0])
b.verify_box()
b.to_csv()