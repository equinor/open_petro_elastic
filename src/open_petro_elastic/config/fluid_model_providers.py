import collections

import pkg_resources

import open_petro_elastic.material.batzle_wang as default
from open_petro_elastic.material.material import vectorize_material


class BatzleWangFluidModelProvider:
    def oil(self, temperature, pressure, reference_density, gas_oil_ratio, gas_gravity):
        return vectorize_material(default.oil)(
            temperature, pressure, reference_density, gas_oil_ratio, gas_gravity
        )

    def gas(self, temperature, pressure, gas_gravity):
        return vectorize_material(default.gas)(temperature, pressure, gas_gravity)

    def brine(self, temperature, pressure, salinity):
        return vectorize_material(default.brine)(temperature, pressure, salinity)

    def condensate(
        self, temperature, pressure, oil_reference_density, gas_oil_ratio, gas_gravity
    ):
        raise NotImplementedError(
            "Default Fluid model provider has no condensate model."
        )


model_names = [
    entry_point.name
    for entry_point in pkg_resources.iter_entry_points(
        "open_petro_elastic.fluid_model_providers"
    )
]

duplicate_names = [
    item for item, count in collections.Counter(model_names).items() if count > 1
]

if duplicate_names:
    raise ImportError(
        "There is a name collision in open_petro_elastic.fluid_model_providers,"
        "check installed plugins.\n"
        f"The following fluid_model_providers have duplicates: ${duplicate_names}"
    )


fluid_model_providers = {
    entry_point.name: entry_point.load()()
    for entry_point in pkg_resources.iter_entry_points(
        "open_petro_elastic.fluid_model_providers"
    )
}

fluid_model_providers["default"] = fluid_model_providers["batzle_wang"]
