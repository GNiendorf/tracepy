## TracePy

[![Build Status](https://travis-ci.org/GNiendorf/tracepy.svg?branch=master)](https://travis-ci.org/GNiendorf/tracepy)
[![Coverage Status](https://coveralls.io/repos/github/GNiendorf/tracepy/badge.svg?branch=master)](https://coveralls.io/github/GNiendorf/tracepy?branch=master)

Ray Tracing and Optical Design in Python

## Overview

TracePy is a sequential ray tracing package written in Python for designing optical systems in the geometric optics regime. It features lens optimization from Scipy and (very soon) Scikit-learn. TracePy is currently in active development and any collaborators would be welcome.

## Installation

To use TracePy I suggest cloning the repository and using the command "pip3 install ." in the downloaded directory. You can also download TracePy directly through pypi with the command below.

```
pip3 install tracepy
```

## Examples

To get started using the software, I suggest looking at the examples provided in the example_projects, optimization_examples, and validation files. The UI for TracePy is most similar to BEAM4, and TracePy's ray tracing algorithm was recreated mostly from Spencer and Murty's iconic paper, "General Ray-Tracing Procedure". Some changes were made to Spencer's algorithm however, and the plotting and optimization software was developed independently by myself.

## Contributing

Contact me at gavinniendorf@gmail.com to get involved in TracePy's development or for questions on how to use it.
