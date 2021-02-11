import numpy


def float_vectorize(f):
    return numpy.vectorize(f, otypes=[float])
