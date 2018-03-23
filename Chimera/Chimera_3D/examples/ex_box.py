from Chimera.Chimera_3D import box

x = 10
y = 10
z = 10
spatial_res = 0.5

b = box.Box(evolution_time=50, conduction=True, settling_mode='stokes terminal', radioactivity=True,
            chemistry=True, verbose=True)
b.build(spatial_res=spatial_res, x=x, y=y, z=z)
b.insert_matrix(material='test_matrix', temperature=2000, depth_range=[0.5, 9.5], conductivity=0.05,
                density=2000, viscosity=0.5)
b.insert_object(material='test_object', temperature=2000, radius=0.01, x=5.2, y=5.4, z=1.0, conductivity=0.05,
                density=2500)
# b.insert_object(material='test_object', temperature=2000, radius=0.01, x=5.2, y=5.4, z=1.0, conductivity=0.05,
#                 density=2500)
# b.insert_object(material='test_object', temperature=2000, radius=0.01, x=3.2, y=9.4, z=1.0, conductivity=0.05,
#                 density=2500)
# b.insert_object(material='test_object', temperature=2000, radius=0.01, x=7.2, y=7.4, z=1.0, conductivity=0.05,
#                 density=2500)
# b.insert_object(material='test_object', temperature=2000, radius=0.01, x=1.2, y=3.4, z=1.0, conductivity=0.05,
#                 density=2500)
# b.insert_object(material='test_object', temperature=2000, radius=0.01, x=3.2, y=1.4, z=1.0, conductivity=0.05,
#                 density=2500)
b.insert_boundary(temperature=3000, depth_range=[0.0, 0.5])
b.insert_boundary(temperature=4000, depth_range=[9.5, 10.0])
b.verify_box()
b.update(animate_model=True)
# b.to_csv()
