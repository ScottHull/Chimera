from Chimera.Chimera_dev import settling_modes


def settling(coord, mode, object, radius, material_properties, spatial_sigfigs, heating=True):

    dimension = len(coord)
    if dimension is 3:
        i, j, k = coord[0], coord[1], coord[2]
        nearest_point = (
            round(i, spatial_sigfigs),
            round(j, spatial_sigfigs),
            round(k, spatial_sigfigs)
        )
    elif dimension is 2:
        i, j = coord[0], coord[1]
        nearest_point = (
            round(i, spatial_sigfigs),
            round(j, spatial_sigfigs)
        )

    if mode is 'stokes_terminal':
        velocity = settling_modes.stokes_terminal(
            material_properties=material_properties,
            matrix_material=None,
            radius=radius,
            object=object
        )
    else:
        pass


