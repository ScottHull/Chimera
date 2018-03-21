from Chimera.Chimera_3D import settling_modes, heat


def settling(coord, mode, object, radius, material_properties, spatial_sigfigs, delta_time, heating=True):

    x, y, z = coord[0], coord[1], coord[2]

    if mode is 'stokes_terminal':
        # currently only works in z-direction--can update to add components
        velocity = settling_modes.stokes_terminal(
            material_properties=material_properties,
            matrix_material=None,
            radius=radius,
            object=object
        )

        x = x + (0.0 * delta_time)
        y = y + (0.0 * delta_time)
        z = z + (velocity * delta_time)

        new_coord = (x, y, z)

        viscous_dissipation = None
        if heating is True:
            viscous_dissipation = heat.viscous_dissipation()

        return new_coord, viscous_dissipation


