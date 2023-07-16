import pathlib
from setuptools import find_packages, setup

HERE = pathlib.Path(__file__).parent

VERSION = '1.0'
PACKAGE_NAME = 'PBHBeta'
URL = 'https://github.com/TadeoDGAguilar'
LICENSE = 'GNU GPL'
DESCRIPTION = 'Library to compute abundances to PBHs in different scenarios on Early Universe'
LONG_DESCRIPTION = (HERE / "README.md").read_text(encoding='utf-8')
LONG_DESC_TYPE = "text/markdown"
INSTALL_REQUIRES = [
      'matplotlib',
      'numpy',
      'scipy',
      ]

setup(
    name=PACKAGE_NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type=LONG_DESC_TYPE,
    author=['Tadeo D. GA', 'Luis E. PA'],
    maintainer = 'Tadeo D. Gómez Aguilar',
    maintainer_email = 'tadeo.daguilar@gmail.com',
    url=URL,
    install_requires=INSTALL_REQUIRES,
    license=LICENSE,
    packages=find_packages(),
    include_package_data=True
)
