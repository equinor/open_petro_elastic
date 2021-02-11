# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import GenericRepr, Snapshot


snapshots = Snapshot()

snapshots["test_fluid_clip_saturations_with_vectors 1"] = GenericRepr(
    "array([[1., 0., 0.],\n       [0., 0., 1.]])"
)

snapshots["test_reservoir_fluids 1"] = (
    GenericRepr(
        "Material(shear_modulus=0, bulk_modulus=[2.7897615e+09], primary_velocity=[1655.40766493], secondary_velocity=[0.], density=[1018.02197492])"
    ),
    GenericRepr(
        "Material(shear_modulus=0, bulk_modulus=[1.4443302e+09], primary_velocity=[1327.76571154], secondary_velocity=[0.], density=[819.2634758])"
    ),
    GenericRepr(
        "Material(shear_modulus=[0.], bulk_modulus=[1.4482208e+09], primary_velocity=[1329.55281424], secondary_velocity=[0.], density=[819.2634758])"
    ),
)

snapshots["test_reservoir_fluids_no_gascond 1"] = (
    GenericRepr(
        "Material(shear_modulus=0, bulk_modulus=[2.7897615e+09], primary_velocity=[1655.40766493], secondary_velocity=[0.], density=[1018.02197492])"
    ),
    GenericRepr(
        "Material(shear_modulus=0, bulk_modulus=[1.4443302e+09], primary_velocity=[1327.76571154], secondary_velocity=[0.], density=[819.2634758])"
    ),
    GenericRepr(
        "Material(shear_modulus=0, bulk_modulus=[74210962.66759408], primary_velocity=[581.97686208], secondary_velocity=[0.], density=[219.10718952])"
    ),
)
