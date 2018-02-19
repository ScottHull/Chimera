from Chimera import linear

z = 100.0
spatial_res = 0.5

line = linear.Line(evolution_time=10000)
line.build(spatial_res=spatial_res, z=z)
line.insert_matrix(material='test_material', initial_temp=2000, conductivity=0.05, temp_grad=0.0, depth_range=[1.0, 99.0])
line.insert_boundary(temperature=2020, depth_range=[0.0, 0.5])
line.insert_boundary(temperature=2020, depth_range=[99.5, 100.0])
line.verify_box()
line.update(timestep=False)
# line.to_csv()
