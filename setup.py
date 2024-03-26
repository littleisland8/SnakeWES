from setuptools import setup,find_packages

setup(name='SnakeWES',
  version=1.0,
  description='Snakemake pipeline for WES analysis',
  url='https://github.com/littleisland8/SnakeWES.git',
  requires=['python (>= 3.6)'],
  author='Simone Romagnoli',
  author_email='simone.romagnoli@unifi.it',
  license='LICENSE.txt',
  #dependency_links=['https://github.com/rrwick/Badread/tarball/master#egg=Badread'],
  #install_requires=['pyfaidx', 'pysam', 'pywgsim', 'pybedtools', 'mappy', 'plotly==3.10.0', 'numpy'],
  zip_safe=False,
  packages=find_packages(),
  include_package_data=True,
  entry_points={'console_scripts': ['SnakeWES=SnakeWES.SnakeWES:main']}          
)