from Chimera import box

x = 5
y = 5
z = 5
spatial_res = 0.2

b = box.Box(evolution_time=20, conduction=True, settling_mode='stokes terminal', radioactivity=True,
            chemistry=True, verbose=True)
b.build(spatial_res=spatial_res, x=x, y=y, z=z)
b.insert_matrix(material='test_matrix', temperature=2000, depth_range=[0.2, 4.8])
b.insert_object(material='test_object', temperature=2000, x=1, y=1, z=1)
b.insert_boundary(temperature=3000, depth_range=[0, 0.2])
b.insert_boundary(temperature=4000, depth_range=[4.8, 5])
b.verify_box()
b.to_csv()
