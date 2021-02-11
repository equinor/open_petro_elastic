# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import GenericRepr, Snapshot


snapshots = Snapshot()

snapshots["test_fluids_vectorized 1"] = GenericRepr(
    "Material(shear_modulus=0, bulk_modulus=[1.64734606e+08 1.83184882e+08], primary_velocity=[1144.72360384 1021.36998928], secondary_velocity=[0. 0.], density=[125.71397711 175.59956795])"
)
