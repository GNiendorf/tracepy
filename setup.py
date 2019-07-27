from setuptools import setup

setup(name='codevp',
      version='0.1',
      description='Optical design software for python.',
      url='http://github.com/GNiendorf/codevp',
      author='Gavin Niendorf',
      author_email='gavinniendorf@gmail.com',
      license='MIT',
      packages=['codevp'],
      install_requires=[
          'numpy',
          'matplotlib',
          'pandas',
          'scipy',
          'scikit-learn',
      ],
      )
