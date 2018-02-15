from Chimera import linear

z = 1*10**4
spatial_res = 1.0

line = linear.Line(evolution_time=20)
line.build(spatial_res=spatial_res, z=z)
line.insert_matrix(material='test_material', initial_temp=2000, temp_grad=0.2, depth_range=[2.0, 998.0])
line.insert_boundary(temperature=3000, depth_range=[0.0, 1.0])
line.insert_boundary(temperature=3000, depth_range=[999.0, 1000.0])
line.verify_box()
line.to_csv()