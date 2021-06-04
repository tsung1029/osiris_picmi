#!/usr/bin/env python

# setup script for osiris_picmi


"""
setup.py file for OSIRIS_PICMI
"""

import argparse
import os
import sys

from setuptools import setup



# Allow to control options via environment vars.
# Work-around for https://github.com/pypa/setuptools/issues/1712
OSPICMI_LIB_DIR = os.environ.get('OSPICMI_LIB_DIR')

setup(name = 'osiris_picmi',
      version = '1.0',
      packages = ['osiris_picmi'],
      package_dir = {'osiris_picmi': 'osiris_picmi'},
      author = 'Frank S. Tsung, Ben J. Winjum, et al'
      author_email ='tsung@physics.ucla.edu'
      description = """PICMI Front End for OSIRIS 4.0""",
      package_data = package_data,
      install_requires = ['numpy', 'picmistandard==0.0.14', 'periodictable'],
      python_requires = '>=3.6',
      zip_safe=False
)