import pandas as pd
from pydantic import root_validator
from pydantic.dataclasses import dataclass

from .dry_rock import DryRock
from .fluids import Fluids
from .minerals import Minerals
from .pressure import Pressure
from .pydantic_config import PetroElasticConfig


class ColumnInsertionError(IndexError):
    pass


def insert_from_dataframe(dataframe, config):
    """
    :param config: A config file that have entries of the form
        {"column": COLUMN_NAME} referring to a column in the dataframe.
    :param dataframe: A dataframe with columns that have values referred to in the
        config dictionary.
    :returns: a new config dictionary with entries of the form
        {"column": COLUMN_NAME} replaced with dataframe[COLUMN_NAME].
    """
    if not isinstance(config, dict) and not isinstance(config, list):
        return config
    if isinstance(config, list):
        return [insert_from_dataframe(dataframe, e) for e in config]
    if list(config.keys()) == ["column"]:
        try:
            return dataframe[config["column"]].to_numpy()
        except (KeyError, ValueError, IndexError) as e:
            raise ColumnInsertionError(
                f"Could not insert column {config['column']}: {e}"
            ) from e
    c = config.copy()
    for key, value in c.items():
        c[key] = insert_from_dataframe(dataframe, value)
    return c


@dataclass(config=PetroElasticConfig)
class Input:
    data: pd.DataFrame = None
    minerals: Minerals = Minerals()
    fluids: Fluids = Fluids()
    dry_rock: DryRock = DryRock()
    pressure: Pressure = Pressure()

    @root_validator(pre=True)
    def insert_dataframe(cls, values):
        if "data" in values and values["data"] is not None:
            return insert_from_dataframe(values["data"], values)
        else:
            cls.check_no_column_reference(values)
        return values

    @classmethod
    def check_no_column_reference(cls, values):
        if isinstance(values, list):
            for v in values:
                cls.check_no_column_reference(v)
        if not isinstance(values, dict):
            return
        if list(values.keys()) == ["column"]:
            raise ColumnInsertionError(
                f"config contains column {values['column']} but no data is given"
            )
        else:
            for v in values.values():
                cls.check_no_column_reference(v)
