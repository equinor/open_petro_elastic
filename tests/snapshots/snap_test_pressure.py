# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import GenericRepr, Snapshot


snapshots = Snapshot()

snapshots["test_pressure_dependency_no_pressure_exp_fit 1"] = GenericRepr(
    "Material(shear_modulus=1.0, bulk_modulus=0.9999999999999998, primary_velocity=1.5275252316519465, secondary_velocity=1.0, density=1.0)"
)

snapshots["test_pressure_dependency_no_pressure_poly_fit 1"] = GenericRepr(
    "Material(shear_modulus=1.0, bulk_modulus=0.9999999999999998, primary_velocity=1.5275252316519465, secondary_velocity=1.0, density=1.0)"
)

snapshots["test_pressure_dependency_no_pressure_power_fit 1"] = GenericRepr(
    "Material(shear_modulus=1.0, bulk_modulus=0.9999999999999998, primary_velocity=1.5275252316519465, secondary_velocity=1.0, density=1.0)"
)
