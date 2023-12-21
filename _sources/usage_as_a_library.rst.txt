Usage as a Library
##################

Generally, the :code:`open_petro_elastic.material` module
is used for interfacing with open_petro_elastic as
a python library:

.. doctest::

    >>> from open_petro_elastic.material import Material
    >>> from open_petro_elastic.material.sandstone import hertz_mindlin
    >>>
    >>> mineral = Material(bulk_modulus=1e9, shear_modulus=1e9, density=1000)
    >>> sand = hertz_mindlin(
    ...     mineral,
    ...     porosity=0.4,
    ...     pressure=1e6
    ... )
    >>> print(sand.density)
    600.0

.. autosummary::
   :toctree: _autosummary
   :recursive:

   open_petro_elastic.material
