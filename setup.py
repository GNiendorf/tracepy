from setuptools import setup

setup(name='raypy',
      version='0.1',
      description='Optical design software for python.',
      url='http://github.com/GNiendorf/raypy',
      author='Gavin Niendorf',
      author_email='gavinniendorf@gmail.com',
      license='MIT',
      packages=['raypy'],
      install_requires=[
          'numpy',
          'matplotlib',
          'pandas',
          'scipy',
          'scikit-learn',
      ],
      )
