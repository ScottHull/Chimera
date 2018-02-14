from random import randint
from . import console
import time


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

    if object_type is 'matrix':
        object_id = random_gen(object_identifier='B', id_val=id_val)
        # while object_id in object_ids:
        #     object_id = random_gen(object_identifier='B')
        return object_id
    elif object_type is 'object':
        object_id = random_gen(object_identifier='A', id_val=id_val)
        return object_id
    elif object_type is 'boundary':
        object_id = random_gen(object_identifier='C', id_val=id_val)
        return object_id
    else:
        object_id = random_gen(object_identifier='Z', id_val=id_val)
        return object_id



def predict_index(coord, max_x, max_y, spatial_res, max_z=None, verbose=True):
    i, j, k = coord[0], coord[1], coord[2]
    s = spatial_res
    x, y, z = max_x, max_y, max_z
    if z is not None:
        index = ((i*y*z)/(s**3)) + (((i*y) + (i*z) + (j*z))/(s**2)) + ((i+j)/s) + (k/s)
        index = int(round((((i*y*z)/(s**3)) + (((i*y) + (i*z) + (j*z))/(s**2)) + ((i+j)/s) + (k/s))))
        return index
    else:
        console.error("2D point prediction not implemented!", verbose=verbose)
