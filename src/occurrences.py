"""
This module provides methods for organizing required raster tiles given a set of species occurrence records
"""

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

class Occurrences(GeoDataFrame):
    """
    A GeoDataFrame with additional methods for interacting with environmental rasters
    """

    def __init__(self, gdf_obj: GeoDataFrame):
        super().__init__(gdf_obj)
        self._obj: GeoDataFrame = gdf_obj

    def within(
        self,
        raster: xarray.DataArray) -> GeoDataFrame:
        """
        Filter the GeoDataFrame to only occurrence records that are within the bounds of a raster

        Parameters
        ----------
        raster: xarray.DataArray
            A raster with which to filter the occurrence records. 

        Returns
        -------
        :obj:`GeoDataFrame`
            Species occurrence records within the bounding box of the raster.
        """
        return self._obj.intersection(box(*raster.rio.bounds()))

    def any_within(
        self,
        raster: xarray.DataArray) -> bool:
        """
        Check if any species occurrence records are within the bounds of a raster

        Parameters
        ----------
        raster: xarray.DataArray
            A raster with which to filter the occurrence records.

        Returns
        -------
        bool:
            True if any species occurrence records are within the raster, otherwise False.
        """
        val = len(self.within(raster).index) > 0
        return val

    # def outliers(
    #     self, 
    #     threshold):
    #     """
    #     Identify outliers based on proximity
    #     """
    #     return None

    def thin(
        self,
        resolution: Optional[float] = 0.5,
        thin_proportion: Optional[float] = 0.75,
        minimum: Optional[int] = None) -> GeoDataFrame:
        """
        Remove some species occurrence records from areas with a high density of species occurrences. 
        Works by projecting a grid with resolution `resolution` onto the occurrence records, then removing a 
        proportion (equal to `thin_proportion`) of occurrences in each grid cell. The occurrences to remove are randomly chosen.

        Parameters
        -----------
        resolution: float
            The area over which to calculate how many species occurrence records must be removed to meet the `thin_proportion`.
            This equals the size of each grid cell in which occurrence records are independently thinned. Defaults to 0.5.
        thin_proportion: float
            The proportion of species occurrence records to remove within each grid cell. Defaults to 0.75.
        minimum:
            The minimum number of species occurrence records that should remain in each grid cell. If `thin_proportion` 
            yields fewer species occurrence records than `minimum`, thinning will stop when the minimum number
            of occurrence records is reached. If a cell has fewer records than `minimum`, none of will be removed. Defaults to None.

        Returns
        -------
        :obj:`GeoDataFrame`:
            A thinned set of the original species occurrence records.

        """
        # Create empty grid to overlay on species occurrence records
        grid = make_grid(gdf_template = self._obj, resolution = resolution)

        # Find which coordinates intersect with the occurrence records, then loop through only those in the for-loop


        # print(grid[0:2, 0:2])
        # For each grid cell, randomly select (1-thin_proportion) of them
        new_points_list = []
        for x in range(len(grid.coords['x'])-1):
            for y in range(len(grid.coords['y'])-1):

                print(x, y)
                # print(grid.coords['x'], 'grid coords x: ', grid.coords['x'][x])
                bbox = grid.coords['x'][x].values, grid.coords['y'][y].values, \
                    grid.coords['x'][x+1].values, grid.coords['y'][y+1].values
                bbox = [int(p) for p in bbox]

                points = self.within(grid[x:x+2, y:y+2])
                points = points[~points.geometry.is_empty]
                nrows = points.shape[0]

                if nrows == 0:
                    continue

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


                print(fr'There are {nrows} records within {bbox}. Removed {nrows - new_points.shape[0]}, yielding {new_points.shape[0]}.')
                
                new_points_list.append(new_points)


        ## Create empty grid to project onto the species occurrence records
        # # Get bounding box of the input GeoDataFrame
        # print(points[~points.geometry.is_empty])
        # return grid
        # return self.within(grid)
        return pd.concat(new_points_list, ignore_index=True)

    def rasterize(
        self, 
        bounds: tuple,
        resolution: int,
        threshold: Optional[int] = 1,
        type: Optional[str] = 'Presence') -> xarray.DataArray:
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

        # Create a raster where values = 1 if there is an occurrence in that tile, otherwise value = 0
        xarray.DataArray

    def read_rasters(
        self,
        distance: int,
        rasters: xarray.DataArray) -> xarray.DataArray:

        # Read in raster tiles that are within `distance` of a species occurrence record.
        
        return self

    

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



