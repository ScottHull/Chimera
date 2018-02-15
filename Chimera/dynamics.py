from . import settling_modes



def settling(coord, mode, material_properties, spatial_sigfigs, heating=True):

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
        pass
