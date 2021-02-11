import numpy as np
import pandas as pd
import pytest

from open_petro_elastic.config.input import insert_from_dataframe


def test_insert_from_dataframe_empty():
    assert insert_from_dataframe(pd.DataFrame(), {}) == {}
    assert insert_from_dataframe(pd.DataFrame(), []) == []
    assert insert_from_dataframe(pd.DataFrame(), 3.0) == 3.0


def test_insert_from_dataframe_inserts_column():
    assert np.all(
        insert_from_dataframe(
            pd.DataFrame({"column_name": [1, 2, 3]}), {"column": "column_name"}
        )
        == [1, 2, 3]
    )


def test_missing_raises():
    with pytest.raises(IndexError, match="insert column"):
        insert_from_dataframe(pd.DataFrame(), {"column": "missing column"})
