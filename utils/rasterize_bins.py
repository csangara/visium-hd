import dask.array as da
import numpy as np
from spatialdata import SpatialData
from spatialdata._types import ArrayLike
from spatialdata.models import Image2DModel
from spatialdata.transformations import get_transformation
from spatialdata_io._constants._constants import VisiumHDKeys
from xarray import DataArray

from sparrow.utils._keys import _REGION_KEY,_INSTANCE_KEY
# inspired by https://github.com/scverse/spatialdata/blob/b9eb240aacd0f68133c95b68ba27529e060ca0b6/src/spatialdata/_core/operations/rasterize_bins.py#L28
def _rasterize_bins( sdata: SpatialData, table_layer: str, obs_column: str | list[str] | None)->DataArray:

    obs_column = [ obs_column ] if isinstance(obs_column, str ) else obs_column

    table = sdata[ table_layer ]

    row_key=VisiumHDKeys.ARRAY_ROW
    col_key=VisiumHDKeys.ARRAY_COL

    min_row, min_col = table.obs[row_key].min(), table.obs[col_key].min()
    n_rows, n_cols = table.obs[row_key].max() - min_row + 1, table.obs[col_key].max() - min_col + 1
    y = (table.obs[row_key] - min_row).values
    x = (table.obs[col_key] - min_col).values


    labels_layer=table.obs[ _REGION_KEY ].cat.categories[0]

    if obs_column is None:

        dtype = table.X.dtype

        keys = table.var_names

        shape = (n_rows, n_cols)

        def channel_rasterization(block_id: tuple[int, int, int] | None) -> ArrayLike:

            image: ArrayLike = np.zeros((1, *shape), dtype=dtype)

            if block_id is None:
                return image

            col = table.X[:, block_id[0]]
            bins_indices, data = col.indices, col.data
            image[0, y[bins_indices], x[bins_indices]] = data

            return image

        image = da.map_blocks(
            channel_rasterization,
            chunks=((1,) * len(keys), *shape),
            dtype=np.uint32,
        )

    else:
        keys = obs_column
        image = np.zeros((len(obs_column), n_rows, n_cols))
        image[:, y, x] = table.obs[obs_column].values.T

    element=Image2DModel.parse( data = image, dims = ("c", "y", "x" ), transformations=get_transformation( sdata.labels[ labels_layer  ], get_all=True ), c_coords = keys )

    return element