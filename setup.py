from setuptools import setup

setup(name='pyBatman',
      version='1.0.0',
      description='NMRbox  Active Directory',
      long_description='A pipeline to infer metabolite concentrations from 1D NMR data '
          'using BATMAN (Bayesian AuTomated Metabolite Analyser for NMR)',
      author='joewandy',
      download_url='https://github.com/glasgowcompbio/pyBatman',
      maintainer='NMRbox',
      install_requires=['scipy<=1.0.0', 'numpy<=1.14.0', 'nmrglue<=0.8',
        'pandas','matplotlib','ipython'],
      maintainer_email='support@nmrbox.org',
      packages=[
          'pyBatman'
      ],
      classifiers=['License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2 :: Only']
      )
