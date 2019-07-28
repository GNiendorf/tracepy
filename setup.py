from setuptools import setup

setup(name='tracepy',
      version='0.1',
      description='Optical design software for python.',
      url='http://github.com/GNiendorf/tracepy',
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
