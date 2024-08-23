import collections

import pkg_resources

import open_petro_elastic.material.batzle_wang as default
from open_petro_elastic.material import span_wagner
from open_petro_elastic.material.conversions import celsius_to_kelvin
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

    def carbon_dioxide(self, *args, **kwargs):
        raise NotImplementedError(
            "Default fluid model provider has no carbon dioxide model"
        )


class SpanWagnerFluidModelProvider:
    def oil(self, *args, **kwargs):
        raise NotImplementedError("Span & Wagner fluid model provider has no oil model")

    def gas(self, *args, **kwargs):
        raise NotImplementedError("Span & Wagner fluid model provider has no gas model")

    def brine(self, *args, **kwargs):
        raise NotImplementedError(
            "Span & Wagner fluid model provider has no brine model"
        )

    def condensate(self, *args, **kwargs):
        raise NotImplementedError(
            "Span & Wagner fluid model provider has no condensate model"
        )

    def carbon_dioxide(self, temperature, pressure, interpolate_density):
        temperature = celsius_to_kelvin(temperature)
        pressure /= 1e6
        return span_wagner.carbon_dioxide(
            temperature,
            pressure,
            None,
            force_vapor="auto",
            interpolate=interpolate_density,
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
