from Chimera import linear, plots

z = 1000.0
spatial_res = 2.0

line = linear.Line(evolution_time=10000)
line.build(spatial_res=spatial_res, z=z)
line.insert_matrix(material='test_material', initial_temp=2000, conductivity=0.5, temp_grad=0.0, depth_range=[2.5, 400.0])
line.insert_matrix(material='test_material_2', initial_temp=2200, conductivity=0.01, temp_grad=0.0, depth_range=[400.5, 600.0])
line.insert_matrix(material='test_material_3', initial_temp=2000, conductivity=0.02, temp_grad=0.0, depth_range=[600.5, 997.5])
line.insert_boundary(temperature=2100, depth_range=[0.0, 2.0])
line.insert_boundary(temperature=2100, depth_range=[998.0, 1000.0])
line.verify_box()
line.update(timestep=False)
plots.temperature_distribution(df=line.get_mesh())
# line.to_csv()
