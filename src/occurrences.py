"""
This module provides methods for organizing required raster tiles given a set of species occurrence records
"""
from __future__ import annotations
from typing import Optional, List, Union

import os
import rioxarray
import rasterio
import xarray
from geopandas import GeoDataFrame
import pandas as pd
from shapely import box, Point
import numpy as np
from utils_occurrences import make_grid




class Occurrences(GeoDataFrame, pd.DataFrame):
    """
    A GeoDataFrame with additional methods for interacting with environmental rasters
    """

    def __init__(self, gdf_obj: GeoDataFrame):
        # super().__init__(gdf_obj)
        super().__init__()
        self._obj: GeoDataFrame = gdf_obj

    # These properties must be overwritten to return class Occurrences rather than GeoDataFrame
    @property
    def _constructor(self):
        return Occurrences

    # @property
    # def _constructor_expanddim(self):
    #     return Occurrences

    @property
    def _constructor_sliced(self):
        return Occurrences

    def within(
        self,
        raster: xarray.DataArray) -> Occurrences:
        """
        Check if geometries from a GeoDataFrame are within the bounds of a DataArray

        Parameters
        ----------
        raster: xarray.DataArray
            A DataArray with which to intersect the geometries

        Returns
        -------
        :obj:`pandas.Series`
            Boolean Series describing whether each row intersects with the raster
        """
        return self._obj.intersects(box(*raster.rio.bounds()))

    def any_within(
        self,
        raster: xarray.DataArray) -> bool:
        """
        Check if any geometries from a GeoDataFrame are within the bounds of a DataArray

        Parameters
        ----------
        raster: xarray.DataArray
            A DataArray with which to intersect the geometries

        Returns
        -------
        bool:
            True if any rows intersect with the raster, otherwise False.
        """
        return bool(self[self.within(raster)].shape[0])

    # def outliers(
    #     self, 
    #     threshold):
    #     """
    #     Identify outliers based on proximity
    #     """
    #     return None

    def thin(
        self,
        resolution: float = 5,
        thin_proportion: Union[Optional[float], Optional[str]] = 'auto',
        minimum: Optional[int] = None) -> Occurrences:
        """
        Remove some species occurrence records from areas with a high density of species occurrences. 
        Works by projecting a grid with resolution `resolution` onto the occurrence records, then removing a 
        proportion (equal to `thin_proportion`) of occurrences in each grid cell. The occurrences to remove are randomly chosen.

        Parameters
        -----------
        resolution: float
            The size of each grid cell in which to downsample records by the `thin_proportion`.
            Defaults to 5 km (25km^2)
        thin_proportion: float
            The proportion of species occurrence records to remove within each grid cell. Defaults to 'auto', which is \
                the the number of records in the cell divided by the resolution^2.
        minimum:
            The minimum number of species occurrence records that should remain in each grid cell. If `thin_proportion` 
            yields fewer species occurrence records than `minimum`, thinning will stop when the minimum number
            of occurrence records is reached. If a cell has fewer records than `minimum`, none of will be removed. Defaults to None.

        Returns
        -------
        :obj:`GeoDataFrame`:
            A thinned set of the original species occurrence records.

        """
        if (type(thin_proportion) == str) and thin_proportion not in ['auto', 'binary']:
            return "`thin_proportion` must be a value between 0 and 1, 'auto', or 'binary'"
        if set(self._obj.geom_type) != {'Point'}:
            return "All geometries must be of geom type 'Point'"

        # Create empty grid to overlay on species occurrence records
        grid = make_grid(gdf_template = self._obj, resolution = resolution)

        # For each grid cell, randomly select (1-thin_proportion) of them
        new_points_list = []
        for x in range(len(grid.coords['x']) - 1):
            for y in range(len(grid.coords['y']) - 1):

                # print(grid.coords['x'], 'grid coords x: ', grid.coords['x'][x])
                bbox = grid.coords['x'][x].values, grid.coords['y'][y].values, \
                    grid.coords['x'][x+1].values, grid.coords['y'][y+1].values
                bbox = [int(p) for p in bbox]

                w = self.within(grid[x:x+2, y:y+2])
                points = self._obj[w]
                
                nrows = points.shape[0]
                if nrows == 0:
                    continue

                if thin_proportion == 'auto':
                    # Must change to vary with resolution and nrows, like 0.75 * (nrows / resolution ** 2) or something
                    thin_proportion = 0.75

                # If the minimum value is None, thin the records in every cell using by the thin_proportion
                if minimum is None:
                    new_points = points.sample(frac = 1 - thin_proportion, replace = False)
                # If the number of records after thinning by thin_proportion is above the minimum, thin the records by thin_proportion
                elif nrows * thin_proportion > minimum:
                    new_points = points.sample(frac = 1 - thin_proportion, replace = False)
                # If the number of records after thinning by thin_proportion is equal to or less than the minimum,
                # thin until the minimum is reached
                elif nrows * thin_proportion <= minimum:
                    # If there are fewer records than the minimum, use the number of records 
                    if nrows < minimum:
                        new_points = points.sample(n = nrows, replace = False)
                    # If there are more records then minimum, use the minimum
                    else: 
                        new_points = points.sample(n = minimum, replace = False)

                # print(fr'There are {nrows} records within {bbox}. Removed {nrows - new_points.shape[0]}, yielding {new_points.shape[0]}.')
                
                new_points_list.append(new_points)

        self._obj = pd.concat(new_points_list, ignore_index=True)

        return self._obj

    def rasterize(
        self, 
        resolution: int,
        threshold: Optional[int] = 1,
        type: Optional[str] = 'presence') -> xarray.DataArray:
        """
        Create a raster to represent the number of species occurrences or the presence-absence of species occurrences. 
        Defaults to a binary raster where 1 represents presence and 0 represents absence.

        Parameters
        -----------
        bounds:
            The bounding box (xmin, ymin, xmax, ymax) of the output raster
        resolution:
            The resolution of each cell in the output raster
        threshold:
            The minimum number of occurrences for a species to be considered present in a raster cell. Defaults to 1.
        type:
            Whether the output raster should represent presence-absence values or a count of the number of species occurrence records.
            Use 'Presence' for presence-absence and 'Count' to count the number of occurrence records.
        
        """
        if type not in ['presence', 'count']:
            return "Type must be equal to 'presence' or 'pount'"

        # Create empty grid to overlay on species occurrence records
        grid = make_grid(gdf_template = self._obj, resolution = resolution)
        
        for x in range(len(grid.coords['x']) - 1):
            for y in range(len(grid.coords['y']) - 1):

                # print(grid.coords['x'], 'grid coords x: ', grid.coords['x'][x])
                bbox = grid.coords['x'][x].values, grid.coords['y'][y].values, \
                    grid.coords['x'][x+1].values, grid.coords['y'][y+1].values
                bbox = [int(p) for p in bbox]

                w = self.within(grid[x:x+2, y:y+2])
                points = self._obj[w]

                if type == 'presence':
                    if not points.empty:
                        grid[x, y] = 1
                elif type == 'count':
                    grid[x, y] = points.shape[0]

        return grid

    def which_rasters(
        self,
        distance: int,
        rasters: Union[List[xarray.DataArray], str]) -> xarray.DataArray:

        rasters_list = []

        # Given a list of raster tiles, find which ones are within x distance of a coord
        buffered_gdf = Occurrences(self._obj.buffer(distance/2, cap_style=3))

        # For each raster tile, check if any of the coordinates are within it
        for raster in rasters:
            if type(raster) == str:
                r = rioxarray.open_rasterio(raster)
            else:
                r = raster

            w = buffered_gdf[buffered_gdf.within(r)]

            # dict = {
            #     "raster": raster,
            #     "records": buffered_gdf[w].index.values.tolist()
            # }

            # rasters_list.append(dict)
            if not w.empty:
                rasters_list.append(raster)

            try:
                del r
            except:
                continue

        return rasters_list

    def load_rasters(
        self,
        distance: int,
        rasters: Union[List[xarray.DataArray], str],
        chunks: str = 'auto') -> xarray.DataArray:
        """
        Load all rasters within `distance` of any species occurrence record.
        """
        to_load = self.which_rasters(distance=distance, rasters=rasters)

        return [rioxarray.open_rasterio(file, chunks='auto') for file in to_load]


    

class Raster(xarray.DataArray):
    """
    An xarray DataArray with methods to interact with GeoDataFrames
    """

    def __init__(self, xda: xarray.DataArray):
        super().__init__(xda)
        self._obj: xarray.DataArray = xda

    def contains(
        self,
        gdf: GeoDataFrame) -> xarray.DataArray:
        return box(*self._obj.rio.bounds()).contains(gdf.geometry)

    def any_contains(
        self,
        gdf: GeoDataFrame) -> bool:
        return any(self.contains(gdf))

# def to_dict(
#         self,
#         rasters: Union[xarray.DataArray, List[xarray.DataArray]]) -> dict:
#         """
#         Return a dict showing which rasters are needed for each observation
#         """
#         # try:
#         #     Iterator = iter(rasters)
#         # except TypeError:
#         #     'Only one raster provided'
#         raster_intersect = []

#         for raster in rasters:
#             occs = self._obj.filter_within()
#             ids = [x for x in occs.id]

#             intersects = {
#                 # "raster": raster.name,
#                 "raster": 5,
#                 "occurrence_ids": ids
#             }

#             raster_intersect.append(intersects)
            
#         return raster_intersect



