# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import GenericRepr, Snapshot


snapshots = Snapshot()

snapshots["test_fluid_substitution 1"] = GenericRepr(
    "Material(shear_modulus=[1.], bulk_modulus=[1.94708995], primary_velocity=[1.28070748], secondary_velocity=[0.70710678], density=[2.])"
)
