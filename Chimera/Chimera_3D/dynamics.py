from math import pi, sqrt

def stokes_terminal(density, drag_coeff, density_matrix, viscosity_matrix, radius,
                    current_coord, lower_model, upper_model, gravity):

    x, y, z = current_coord[0], current_coord[1], current_coord[2]
    if upper_model < z < lower_model:  # ensure that the object exists outside of the upper and lower boundary layers
        diameter = 2 * radius
        friction_coeff = (pi/6) * ((density - density_matrix) / density_matrix) * \
                (density_matrix / viscosity_matrix)**2 * (gravity * diameter**3)
        if friction_coeff < 10:  # low frictional coefficient, when body is in laminar flow regime
            velocity = ((density - density_matrix) * gravity * diameter ** 2) / (
                    18 * viscosity_matrix)  # calculates the velocity of the body
            return (0.0, 0.0, velocity)  # x, y, z component velocity
        else:
            velocity = sqrt((((4 * gravity * diameter) / (3 * drag_coeff)) * (((density - density_matrix) / density_matrix))))  # high frictional coefficient, when body is in turbulent flow regime
            return (0.0, 0.0, velocity)  # x, y, z component velocity
    else:
        velocity = (0.0, 0.0, 0.0)
        return velocity

def grav_force(object_radius, object_density, grav_accel):
    """
    The force due to gravity on an object.
    :param mass:
    :param grav_accel:
    :return:
    """
    # fg = mass * grav_accel
    fg = (pi / 6) * ((2 * object_radius)**3) * object_density * grav_accel
    return fg

def buoyant_force(matrix_density, object_radius, grav_accel):
    """
    The buoyant force on an object due to the surrounding matrix.
    :param matrix_density:
    :param object_volume:
    :param grav_accel:
    :return:
    """
    # fb = matrix_density * object_volume * grav_accel
    fb = (pi / 6) * ((2 * object_radius)**3) * matrix_density * grav_accel
    return fb

def drag_force(drag_coeff, matrix_density, object_radius, distance_travelled, delta_time):
    """
    The drag force exerted on an object.
    :param drag_coeff:
    :param matrix_density:
    :param object_velocity:
    :param body_radius:
    :return:
    """
    object_velocity = distance_travelled[2] / delta_time
    fd = drag_coeff * 0.5 * matrix_density * (object_velocity**2) * (pi * (object_radius**2))
    return fd

def work_conservative(force, velocity, distance_travelled, delta_time):
    """
    The work done assuming a conservative force, W = F * d.
    :param force:
    :param velocity:
    :param delta_time:
    :return:
    """
    w = force * distance_travelled[2]
    return w
