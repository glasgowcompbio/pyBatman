#!/usr/bin/env python3.6
from setuptools import setup

setup(name='pyBatman',
      version='1.0.0',
      description='NMRbox  Active Directory',
      maintainer='nmrbox team',
      install_requires=['scipy==1.0.0','numpy==1.14.0','nmrglue'],
      maintainer_email='support@nmrbox.org',
      packages=[
          'pyBatman'
      ]
      )
