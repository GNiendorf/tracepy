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

To get started using the software, I suggest looking at the examples provided in the examples folder. The UI for TracePy is most similar to BEAM4, and TracePy's ray tracing algorithm was recreated mostly from Spencer and Murty's iconic paper, "General Ray-Tracing Procedure". Some changes were made to Spencer's algorithm however, and the plotting and optimization software was developed independently by myself.

## Contributing

I would suggest reading Spencer and Murty's paper, ["General Ray-Tracing Procedure"](https://www.osapublishing.org/josa/viewmedia.cfm?uri=josa-52-6-672&seq=0) for an overview of TracePy's algorithm. There are also several open issues that we welcome collaborators on. If you have any questions feel free to contact me at gavinniendorf@gmail.com.
