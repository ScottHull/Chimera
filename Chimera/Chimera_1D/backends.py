from . import console

def generate_object_id(object_type, id_val):
    """
    Generates object ID codes so that specific objects and materials can be tracked
    Object/matrial types are unique and coded by the a letter specifying the general type followed by a
    unique number combination
    The general object/material types are coded by the first letter as follows:
        - 'A' = object
        - 'B' = matrix
        - 'C' = boundary
    :param matrix:
    :return: object_id
    """
    def random_gen(object_identifier, id_val):
        """
        Generates the ID.  Will be called continuously until a unique ID is returned.
        :param object_identifier:
        :param id_val:
        :return:
        """
        object_id = object_identifier + str(id_val)

        return object_id

    object_type = object_type.lower()

    if object_type == 'matrix':
        object_id = random_gen(object_identifier='B', id_val=id_val)
        return object_id
    elif object_type == 'object':
        object_id = random_gen(object_identifier='A', id_val=id_val)
        return object_id
    elif object_type == 'boundary':
        object_id = random_gen(object_identifier='C', id_val=id_val)
        return object_id
    else:
        object_id = random_gen(object_identifier='Z', id_val=id_val)
        return object_id

def predict_index(coord, max_x, max_y, spatial_res, max_z=None, verbose=True):
    """
    Predicts the index position of a coordinate in three dimensional (3D) space.
    :param coord:
    :param max_x:
    :param max_y:
    :param spatial_res:
    :param max_z:
    :param verbose:
    :return:
    """
    i, j, k = coord[0], coord[1], coord[2]
    s = spatial_res
    x, y, z = max_x, max_y, max_z
    if z is not None:
        index = int(round((((i*y*z)/(s**3)) + (((i*y) + (i*z) + (j*z))/(s**2)) + ((i+j)/s) + (k/s))))
        return index
    else:
        console.error("2D point prediction not implemented!", verbose=verbose)

def predict_linear_index(z, spatial_res):
    """
    Returns the index position of a coordinate specified in a one dimensional (1D) model.
    :param z:
    :param spatial_res:
    :return:
    """
    index = int(round(z / spatial_res))
    return index

def override_timestep(timestep, conductivities, spatial_res, spatial_sigfigs):
    """
    Sets the time intervals for model evolution.  An override allows the user to define a custom time interval
    with the risk of poor convergence.
    :param override:
    :return:
    """
    if timestep is False:
        # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        # the explicit method timestep is only stable if:
        # deltaT <= (deltaX^2)/(2*(conductivity))
        # stability requires the maximum usage of conductivity to get the appropriate timestep
        # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        max_conductivity = max(conductivities)
        delta_time = (spatial_res**2) / (2 * (max_conductivity))
        delta_time = round(delta_time, spatial_sigfigs)
        return delta_time
    else:
        delta_time = timestep
        return delta_time

def calculate_thermal_diffusivity(thermal_conductivity, density, specific_heat_capacity):
    """
    Returns thermal diffusivity from thermal conductivity.
    alpha = k/(rho*c_p)
    :param thermal_conductivity:
    :param density:
    :param specific_heat_capacity:
    :return:
    """
    diffusivity = thermal_conductivity / (density * specific_heat_capacity)
    return diffusivity

def calculate_thermal_conductivity(thermal_diffusivity, density, specific_heat_capacity):
    """
    Returns thermal conductivity from thermal diffusivity.
    k = alpha * c_p * rho
    :param thermal_diffusivity:
    :param density:
    :param specific_heat_capacity:
    :return:
    """
    conductivity = thermal_diffusivity * density * specific_heat_capacity
    return conductivity

def convect_heat_transfer_coeff(heat_flux, temp_top, temp_bottom):
    """
    The heat transfer coefficient (or film coefficient/efficiency) is the ratio of heat flux to the temperature
    gradient between the solid surface and the surrounding fluid area, k.
    h = q/(delta_T), where:
    q = heat flux
    delta_T = temperature difference between the solid surface and the surrounding fluid area, k
    :param heat_flux:
    :param temp_difference:
    :return:
    """
    q = heat_flux / (temp_top - temp_bottom)
    return q

def gravitational_acceleration(radius, mass):
    pass
