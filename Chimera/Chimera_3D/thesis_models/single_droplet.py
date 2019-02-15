from Chimera.Chimera_3D import box
import random

# define some spatial parameters used to build the box model
x = 1200.0
y = 1200.0
z = 262.0 * 1000
spatial_res = 400.0

# define some material parameters used to build the box model
density_silicate = 3750
cp_silicate = 10**3
thermal_diffusivity_silicate = 10**-6
conductivity_silicate = thermal_diffusivity_silicate * density_silicate * cp_silicate
# conductivity_silicate = 0.2
viscosity_silicate = 10**-3.5

density_liq_fe = 7800
cp_liq_fe = 825
conductivity_liq_fe = 0.45
drag_coeff_liq_fe = 0.2

w_chemical_diffusivity = 10**-8

# instantiate the box
b = box.Box(
    evolution_time=1015200,
    gravity=0.25,
    conduction=False,
    settling_mode='stokes terminal',
    radioactivity=True,
    chem=True,
    verbose=True,
    multiprocessing=False,
    num_processors=1,
)
b.build(
    spatial_res=spatial_res,
    x=x,
    y=y,
    z=z
)
b.insert_matrix(
    material='silicate_melt',
    depth_range=[spatial_res, z - 400],
    conductivity=conductivity_silicate,
    density=density_silicate,
    viscosity=viscosity_silicate,
    heat_capacity=cp_silicate,
    composition={'w': 100},
    element_diffusivity={'w': w_chemical_diffusivity},
    temperature=2000,
    pressure=0,
    fO2=-2.45,
    grad_temperature=0.012,
    grad_pressure=(3.75 * 10**(-4)),
    grad_fO2=0,
)
for i in range(1):
    b.insert_object(
        material='test_object',
        temperature=2000,
        radius=0.0185,
        x=600,
        y=600,
        z=401,
        conductivity=conductivity_liq_fe,
        density=density_liq_fe,
        drag_coeff=drag_coeff_liq_fe,
        cp=cp_liq_fe,
        composition={'w': 10**-6}
    )

b.insert_boundary(
    temperature=2000,
    depth_range=[0.0, spatial_res],
    location='top',
    composition={'w': 100.0},
    element_diffusivity={'w': w_chemical_diffusivity}
)
b.insert_boundary(
    temperature=2008,
    depth_range=[z - spatial_res, z],
    location='bottom',
    composition={'w': 100.0},
    element_diffusivity={'w': w_chemical_diffusivity}
)

b.verify_box()
b.update(
    animate_model=False,
    show_model=False,
    timestep=1549.97
)

b.to_csv()