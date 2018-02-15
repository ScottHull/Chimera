from math import pi, sqrt


def stokes_terminal(material_properties, radius, matrix_material, object):

    diameter = radius * 2
    density_body = material_properties['density'][object]
    density_matrix = material_properties['density'][matrix_material]
    drag_coeff = material_properties['drag coefficient'][object]
    matrix_viscosity = material_properties['viscosity'][matrix_material]
    gravity = 9.81
    friction_coeff = (pi/6) * ((density_body - density_matrix) / density_matrix) * \
            (density_matrix / matrix_viscosity)**2 * (gravity * diameter**3)

    if friction_coeff < 10:  # low frictional coefficient, when body is in laminar flow regime
        velocity = ((density_body - density_matrix) * gravity * diameter ** 2) / (
                18 * matrix_viscosity)  # calculates the velocity of the body
        return velocity
    else:
        velocity = sqrt(((4 / (3 * drag_coeff)) * (((density_body - density_matrix) / density_matrix) * (gravity * diameter))))
        return velocity
