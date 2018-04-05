from math import pi, sqrt
from Chimera.Chimera_3D import backends


# def stokes_terminal(density, drag_coeff, density_matrix, viscosity_matrix, radius,
#                     current_coord, lower_model, upper_model):
#
#     x, y, z = current_coord[0], current_coord[1], current_coord[2]
#     if upper_model < z < lower_model:  # ensure that the object exists outside of the upper and lower boundary layers
#         diameter = radius * 2
#         gravity = 9.81
#         friction_coeff = (pi/6) * ((density - density_matrix) / density_matrix) * \
#                 (density_matrix / viscosity_matrix)**2 * (gravity * diameter**3)
#
#         if friction_coeff < 10:  # low frictional coefficient, when body is in laminar flow regime
#             velocity = ((density - density_matrix) * gravity * diameter ** 2) / (
#                     18 * viscosity_matrix)  # calculates the velocity of the body
#             return (0.0, 0.0, velocity)  # x, y, z component velocity
#         else:
#             velocity = sqrt(((4 / (3 * drag_coeff)) * (((density - density_matrix) / density_matrix) *
#                             (gravity * diameter))))  # high frictional coefficient, when body is in turbulent flow regime
#             return (0.0, 0.0, velocity)  # x, y, z component velocity
#     else:
#         velocity = (0.0, 0.0, 0.0)
#         return velocity
