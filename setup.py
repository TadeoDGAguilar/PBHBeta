import pathlib
from setuptools import find_packages, setup

HERE = pathlib.Path(__file__).parent

VERSION = '1.0.5'
PACKAGE_NAME = 'PBHBeta'
URL = 'https://github.com/TadeoDGAguilar/PBHBeta'
LICENSE = 'GNU GPL'
DESCRIPTION = 'PBHBeta, a code to explore the constraints on PBH abundances and power spectrum.'
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
