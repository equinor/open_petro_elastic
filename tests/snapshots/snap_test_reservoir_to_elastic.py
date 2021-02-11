# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import GenericRepr, Snapshot


snapshots = Snapshot()

snapshots["test_reservoir_to_elastic 1"] = (
    {
        "aisat": GenericRepr("array([4359286.49764948])"),
        "ksat": GenericRepr("array([5.52268254e+09])"),
        "mysat": GenericRepr("array([3.38138982e+09])"),
        "rsat": GenericRepr("array([1894.42682868])"),
        "sisat": GenericRepr("array([2530967.32511557])"),
        "vpsat": GenericRepr("array([2301.11104407])"),
        "vpvssat": GenericRepr("array([1.7223796])"),
        "vssat": GenericRepr("array([1336.00690552])"),
    },
    {
        "kdry": GenericRepr("array([5.52268254e+09])"),
        "mydry": GenericRepr("array([3.38138982e+09])"),
        "rdry": GenericRepr("array([1696.])"),
        "vpdry": GenericRepr("array([2432.])"),
        "vsdry": GenericRepr("array([1412.])"),
    },
    {
        "kdryref": GenericRepr("array([5.52268254e+09])"),
        "mydryref": GenericRepr("array([3.38138982e+09])"),
        "rdryref": GenericRepr("array([1696.])"),
        "vpdryref": GenericRepr("array([2432.])"),
        "vsdryref": GenericRepr("array([1412.])"),
    },
    {
        "kdryphi": 5522682538.666667,
        "mydryphi": 3381389824.0,
        "rdryphi": 1696.0,
        "vpdryphi": 2432.0,
        "vsdryphi": 1412.0,
    },
    {
        "kmin": 36.8,
        "kmin_fls": GenericRepr("array([5.52268254e+09])"),
        "mymin": 44.0,
        "rmin": 2650.0,
    },
    {
        "kfl": GenericRepr("array([4.99107669e+09])"),
        "rfl": GenericRepr("array([551.18563522])"),
    },
    {
        "kgas": GenericRepr("array([1.4482208e+09])"),
        "koil": GenericRepr("array([1.4443302e+09])"),
        "kwat": GenericRepr("array([2.7897615e+09])"),
        "rgas": GenericRepr("array([819.2634758])"),
        "roil": GenericRepr("array([819.2634758])"),
        "rwat": GenericRepr("array([1018.02197492])"),
    },
)
