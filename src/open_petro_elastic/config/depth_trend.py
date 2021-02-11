from typing import Optional

import numpy as np
from pydantic.dataclasses import dataclass
from typing_extensions import Literal

from .coefficients import PolyfitCoefficients
from .pydantic_config import PetroElasticConfig
from .vector_type import Array


@dataclass(config=PetroElasticConfig)
class DepthTrend:
    """
    Adjusts given dry rock to depth. Difference is
    calculated from given polynomial.
    :param reference_dry_rock: Dry rock to adjust at reference depth.
    :param depth: Difference in depth from reference.
    :param coefficients: The polynomial coefficients to used to adjust the material,
    such that the terms in the polynomial are p[i,j] * depth**i *
    mineral.bulk_modulus ** j.
    """

    depth: Array[float]
    coefficients: PolyfitCoefficients
    type: Literal["depth_trend"] = "depth_trend"
    reference_depth: Optional[Array[float]] = None
    max_depth: Optional[Array[float]] = None

    @property
    def effective_depth(self):
        result = np.asarray(self.depth).copy()
        if self.max_depth is not None:
            np.clip(result, None, self.max_depth)
        if self.reference_depth is not None:
            result = result - self.reference_depth
        return result

    def adjust_material(
        self,
        material,
        pressure=None,
        porosity=None,
        mineral=None,
    ):
        """
        :returns: The material at adjusted depth.
        """
        return self.coefficients.polyfit(
            self.effective_depth,
            material,
        )
