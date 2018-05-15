from Chimera.Chimera_3D import box
import random

x = 5
y = 5
z = 15
spatial_res = 0.5

density_silicate = 3750
cp_silicate = 10**3
diffusivity_silicate = 10**-6
conductivity_silicate = diffusivity_silicate * density_silicate * cp_silicate
# conductivity_silicate = 0.2
viscosity_silicate = 10**-2

density_liq_fe = 7800
cp_liq_fe = 825
conductivity_liq_fe = 0.45
drag_coeff_liq_fe = 0.2

b = box.Box(
    evolution_time=50,
    conduction=True,
    settling_mode='stokes terminal',
    radioactivity=True,
    chem=True,
    verbose=True,
    multiprocessing=False,
    num_processors=1
)
b.build(
    spatial_res=spatial_res,
    x=x,
    y=y,
    z=z
)
b.insert_matrix(
    material='test_matrix',
    depth_range=[round(0.0 + spatial_res, 2), round(z - spatial_res, 2)],
    conductivity=conductivity_silicate,
    density=density_silicate,
    viscosity=viscosity_silicate,
    heat_capacity=cp_silicate,
    composition={'w': 100},
    element_diffusivity={'w': 1.5 * 10**-2},
    temperature=2000,
    pressure=1,
    fO2=-1.0,
    grad_temperature=0.5,
    grad_pressure=0,
    grad_fO2=-0.001
)
for i in range(1):
    b.insert_object(
        material='test_object',
        temperature=2000,
        radius=0.01,
        x=random.uniform(0, x),
        y=random.uniform(0, y),
        z=random.uniform(spatial_res * 3, spatial_res * 4),
        conductivity=conductivity_liq_fe,
        density=density_liq_fe,
        drag_coeff=drag_coeff_liq_fe,
        cp=cp_liq_fe,
        composition={'w': 1.0}
    )
# b.insertObjectGenerator(generator_name='test',
#                         num_per_iteration=5,
#                         material='test_gen_object',
#                         temperature=2000,
#                         radius=0.01,
#                         x_range=[0, x],
#                         y_range=[0, y],
#                         z_range=[spatial_res * 3, spatial_res * 3 + spatial_res],
#                         conductivity=conductivity_liq_fe,
#                         density=density_liq_fe,
#                         drag_coeff=drag_coeff_liq_fe,
#                         cp=cp_liq_fe,
#                         composition={'w': 1.0})
b.insert_boundary(
    temperature=2000,
    depth_range=[0.0, round(0.0 + spatial_res), 2],
    location='top',
    composition={'w': 100.0},
    element_diffusivity={'w': 1.5 * 10**-2}
)
b.insert_boundary(
    temperature=2000,
    depth_range=[round(z - spatial_res, 2), z],
    location='bottom',
    composition={'w': 100.0},
    element_diffusivity={'w': 1.5 * 10**-2}
)
b.verify_box()
b.update(
    animate_model=False,
    show_model=False,
    timestep=1.0
)
# b.to_csv()
