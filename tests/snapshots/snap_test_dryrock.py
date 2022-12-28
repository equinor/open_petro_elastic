# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import GenericRepr, Snapshot


snapshots = Snapshot()

snapshots["test_dryrock_phi_fitted_vp_vs 1"] = (
    GenericRepr(
        "Material(shear_modulus=[7.6585e+09], bulk_modulus=[1.20751667e+10], primary_velocity=[2900.], secondary_velocity=[1700.], density=[2650.])"
    ),
    GenericRepr(
        "Material(shear_modulus=44.0, bulk_modulus=36.8, primary_velocity=0.18980294316133353, secondary_velocity=0.128855630784633, density=2650.0)"
    ),
)

snapshots["test_dryrock_phi_k_my 1"] = (
    GenericRepr(
        "Material(shear_modulus=[12.], bulk_modulus=[20.], primary_velocity=[0.1165543], secondary_velocity=[0.06729266], density=[2650.])"
    ),
    GenericRepr(
        "Material(shear_modulus=44.0, bulk_modulus=36.8, primary_velocity=0.18980294316133353, secondary_velocity=0.128855630784633, density=2650.0)"
    ),
)
