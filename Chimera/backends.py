from Chimera.Chimera_dev import console


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
        object_id = object_identifier + str(id_val)

        return object_id

    object_type = object_type.lower()

    if object_type == 'matrix':
        object_id = random_gen(object_identifier='B', id_val=id_val)
        # while object_id in object_ids:
        #     object_id = random_gen(object_identifier='B')
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
        # the explicit method timestep is only stable if:
        # deltaT <= (deltaX^2)/(2*(conductivity))
        # stability requires the maximum usage of conductivity to get the appropriate timestep
        max_conductivity = max(conductivities)
        delta_time = (spatial_res**2) / (2 * (max_conductivity))
        delta_time = round(delta_time, spatial_sigfigs)
        return delta_time
    else:
        delta_time = timestep
        return delta_time


