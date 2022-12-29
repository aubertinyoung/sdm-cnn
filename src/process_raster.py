"""
This module provides functions for organizing rasters
"""

from typing import Optional

import os
import rioxarray
import rasterio
import xarray
from geopandas import GeoDataFrame
from shapely import box

def chunk_raster(
    raster: xarray.DataArray,
    resolution: int,
    bbox: Optional[tuple] = None,
    output_path: Optional[str] = os.path.abspath(os.curdir)) -> xarray.DataArray:
    """
    Split a large raster file into many smaller files with equal-sized raster tiles

    Parameters
    -----------
    raster: :obj:`xarray.DataArray`
        The large raster file or the path to the large raster file
    res: int
        The approximate resolution of the output tiles. The true output resolution may differ in order to make tiles have equal size.
    output_path: str, optional
        File location in which the output tiles will be written to. Defaults to your current directory.
    """

    if type(raster) == str:
        r = rioxarray.open_rasterio(raster)
    else:
        r = raster

    # Get bounds of the raster, then split into equal-sized tiles along both x and y axes
    if bbox is None:
        xmin, ymin, xmax, ymax = r.rio.bounds()
    elif len(bbox) == 4:
        xmin, ymin, xmax, ymax = bbox
    else:
        return '`bbox` requires four values: xmin, ymin, xmax, ymax.'

    y_length = ymax - ymin # Length of the y axis
    x_length = xmax - xmin # Length of the x axis
    tiles_y = int(y_length // resolution) # Number of even-sized tiles along the y axis that can be made with the given res
    tiles_x = int(x_length // resolution) # Number of even-sized tiles along the x axis that can be made with the given res

    # If the x length is 360, dividing by a res of 2 returns 180. So the length of each tile should be 180
    # print(tiles_x)

    for xdx, x in enumerate(range(tiles_x)):
        for ydx, y in enumerate(range(tiles_y)):
            xmin_t = xmin + xdx * resolution
            xmax_t = xmin_t + resolution
            ymin_t = ymin + ydx * resolution
            ymax_t = ymin_t + resolution
            file_name = f'{round(xmin_t, 2)}_{round(ymin_t, 2)}_{round(xmax_t, 2)}_{round(ymax_t, 2)}'

            try:
                tile = r.rio.clip_box(xmin_t, ymin_t, xmax_t, ymax_t)
                tile.rio.to_raster(f'{output_path}/{file_name}.tif')
                # print(f'Created raster for {file_name} at {output_path}')
                del tile
            except:
                continue
            
    return None

# def load_rasters(
#     rasters,
#     dims
# )

# class Raster(xarray.DataArray):
#     """
#     Creates a new class for xarray.DataArray objects with additional methods
#     """
#     def __init__(self):
#         return self

#     def contains_sp(self, occs: GeoDataFrame) -> bool:
#         """
#         Check if the tile contains any species occurrence records in a GeoDataFrame
#         """
#         return any(occs.intersection(box(*self.rio.bounds())))
