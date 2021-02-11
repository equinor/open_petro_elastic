# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import GenericRepr, Snapshot


snapshots = Snapshot()

snapshots["test_dryrock_depthtrend_k_my 1"] = GenericRepr(
    "Material(shear_modulus=[12.], bulk_modulus=[20.], primary_velocity=[0.1165543], secondary_velocity=[0.06729266], density=[2650.])"
)

snapshots["test_dryrock_depthtrend_vp_vs 1"] = GenericRepr(
    "Material(shear_modulus=[12.], bulk_modulus=[20.], primary_velocity=[0.1165543], secondary_velocity=[0.06729266], density=[2650.])"
)
