from setuptools import setup

setup(name='RayPy',
      version='0.1',
      description='Optical design software for python.',
      url='http://github.com/GNiendorf/RayPy',
      author='Gavin Niendorf',
      author_email='gavinniendorf@gmail.com',
      license='MIT',
      packages=['RayPy'],
      install_requires=[
          'numpy',
          'matplotlib',
          'pandas',
          'scipy',
          'scikit-learn',
      ],
      )
