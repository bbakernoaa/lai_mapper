try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(name='lai_mapper',
      version='0.0.1',
      url='https://github.com/bbakernoaa/lai_mapper',
      license='GNU',
      include_package_data=True,
      author='Barry D. Baker',
      author_email='barry.baker@noaa.gov',
      maintainer='Barry Baker',
      maintainer_email='barry.baker@noaa.gov',
      packages=find_packages(),
      keywords=['lai', 'model', 'satellite'],
      description='LAI global map generator from satellites',
      install_requires=[
          'netcdf4', 'xarray', 'numpy', 'argparse', 'pandas', 'rasterio'
      ])
