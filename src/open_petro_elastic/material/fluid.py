from .material import Material


def fluid(*args, **kwargs):
    """
    fluid is a material with shear modulus=0
    """
    return Material(shear_modulus=0, *args, **kwargs)
