from Chimera import linear, plots

z = 10000.0  # currently in m
spatial_res = 5.0  # currently in m

# instantiate the linear model
line = linear.Line(evolution_time=60*60*24, time_units='seconds')

# construct the linear model
line.build(spatial_res=spatial_res, z=z)

# insert matrix materials
line.insert_matrix(material='test_material_1', initial_temp=2800, conductivity=0.8, temp_grad=0.1, depth_range=[10.0, 4000.0])
line.insert_matrix(material='test_material_2', initial_temp=2000, conductivity=0.01, temp_grad=0.1, depth_range=[4005.0, 6000.0])
line.insert_matrix(material='test_material_3', initial_temp=1000, conductivity=0.8, temp_grad=0.1, depth_range=[6005.0, 7000.0])
line.insert_matrix(material='test_material_4', initial_temp=2000, conductivity=0.01, temp_grad=0.1, depth_range=[7005.0, 8000.0])
line.insert_matrix(material='test_material_5', initial_temp=2800, conductivity=0.8, temp_grad=0.1, depth_range=[8005.0, 9990.0])

# insert boundary layers
line.insert_boundary(temperature=2800, depth_range=[0.0, 5.0])
line.insert_boundary(temperature=2800, depth_range=[9995.0, 10000.0])

# verify that the box as been constructed properly
line.verify_box()

# run the model
# timestep=False parameter automatically calculates the appropriate timestep
# auto_update=False parameter automatically runs the model to the end of the defined evolution time
line.update(timestep=False)

# plot the temperature & nusselt number distribution of the model
plots.temperature_distribution(mesh=line.get_mesh())
plots.nusselt_distrbution(layers=line.get_layers(), mesh=line.get_mesh())

# write the mesh data to a csv file
# line.to_csv()
