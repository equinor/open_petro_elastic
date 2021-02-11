from dataclasses import field

import numpy as np
from pydantic import conlist
from pydantic.dataclasses import dataclass

from open_petro_elastic.float_vectorize import float_vectorize
from open_petro_elastic.material import Material, hashin_shtrikman_average

from .constituent import Constituent, fix_one_fraction
from .pydantic_config import PetroElasticConfig

DEFAULT_SAND = Material(density=2650.0, bulk_modulus=36.8e9, shear_modulus=44.0e9)
DEFAULT_SHALE = Material(density=2650.0, bulk_modulus=25.0e9, shear_modulus=12.0e9)


@dataclass(config=PetroElasticConfig)
class Minerals:
    """
    Minerals is a composite material of several mineral materials.
    :param constituents: The list of constituent minerals.
    """

    constituents: conlist(Constituent, min_items=1) = field(
        default_factory=(
            lambda: [
                Constituent(DEFAULT_SAND, 0.5),
                Constituent(DEFAULT_SHALE, 0.5),
            ]
        )
    )

    def __len__(self):
        return len(self.constituents)

    def __iter__(self):
        return iter(self.constituents)

    def __getitem__(self, item):
        return self.constituents[item]

    def __post_init_post_parse__(self):
        fix_one_fraction(self)
        if not np.all(np.isclose(sum([c.fraction for c in self]), 1.0, atol=1e-5)):
            raise ValueError("Sum of fractions must be one")

    @property
    def as_mixture(self):
        """
        :returns: A material representing the mixture of constituents.
        """
        if len(self) == 0:
            raise ValueError("Cannot mix empty material list")
        mixed_mineral = self[0]

        @float_vectorize
        def part_of_total(fraction, total):
            if total == 0:
                return 1.0
            return fraction / total

        for next_mineral in self[1:]:
            sub_total = mixed_mineral.fraction + next_mineral.fraction
            fraction = part_of_total(mixed_mineral.fraction, sub_total)
            new_mix = hashin_shtrikman_average(
                mixed_mineral.material, next_mineral.material, fraction
            )
            mixed_mineral = Constituent(
                material=new_mix,
                fraction=sub_total,
            )
        return mixed_mineral.material
