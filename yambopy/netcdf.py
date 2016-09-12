#we try to use netcdf, if not present we won't use it
try:
    from netCDF4 import Dataset
    _has_netcdf = True
except ImportError:
    _has_netcdf = False
