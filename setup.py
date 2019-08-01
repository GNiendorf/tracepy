from setuptools import setup

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

setup(name='tracepy',
      version='0.1.1',
      description='Optical design software for python.',
      url='http://github.com/GNiendorf/tracepy',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Gavin Niendorf',
      author_email='gavinniendorf@gmail.com',
      license='MIT',
      packages=['tracepy'],
      install_requires=[
          'numpy',
          'matplotlib',
          'pandas',
          'scipy',
          'scikit-learn',
      ],
      )
