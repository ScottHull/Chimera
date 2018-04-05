from Chimera.Chimera_3D import box
import random

x = 10
y = 10
z = 20
spatial_res = 0.5

b = box.Box(evolution_time=50, conduction=True, settling_mode='stokes terminal', radioactivity=True,
            chemistry=True, verbose=True, multiprocessing=False, num_processors=2)
b.build(spatial_res=spatial_res, x=x, y=y, z=z)
b.insert_matrix(material='test_matrix', temperature=2000, depth_range=[round(0.0 + spatial_res, 2), round(z-spatial_res, 2)], conductivity=0.05,
                density=2000, viscosity=0.5)
for i in range(1):
    b.insert_object(material='test_object', temperature=2000, radius=0.01, x=random.uniform(0, x),
                    y=random.uniform(0, y), z=random.uniform(spatial_res + 0.001, z - (2 * spatial_res)),
                    conductivity=0.05, density=7800, drag_coeff=0.2, cp=0.45)
# b.insert_object(material='test_object', temperature=2000, radius=0.01, x=5.2, y=5.1, z=1.0, conductivity=0.05,
#                 density=2500)
# b.insert_object(material='test_object', temperature=2000, radius=0.01, x=5.2, y=5.9, z=1.0, conductivity=0.05,
#                 density=2500)
# b.insert_object(material='test_object', temperature=2000, radius=0.01, x=3.2, y=9.2, z=1.0, conductivity=0.05,
#                 density=2500)
# b.insert_object(material='test_object', temperature=2000, radius=0.01, x=7.8, y=7.4, z=1.0, conductivity=0.05,
#                 density=2500)
# b.insert_object(material='test_object', temperature=2000, radius=0.01, x=1.2, y=3.8, z=1.0, conductivity=0.05,
#                 density=2500)
# b.insert_object(material='test_object', temperature=2000, radius=0.01, x=3.6, y=1.4, z=1.0, conductivity=0.05,
#                 density=2500)
b.insert_boundary(temperature=2000, depth_range=[0.0, round(0.0 + spatial_res), 2], location='top')
b.insert_boundary(temperature=2000, depth_range=[round(z - spatial_res, 2), z], location='bottom')
b.verify_box()
b.update(animate_model=False)
b.to_csv()
