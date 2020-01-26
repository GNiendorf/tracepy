## <img alt="TracePy" src="https://user-images.githubusercontent.com/25272611/62305283-dc62a300-b43c-11e9-8436-d88c8555b110.png" height="60">

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/454c6504e63f4accaa9353e7dcfda00e)](https://app.codacy.com/app/gavinniendorf/tracepy?utm_source=github.com&utm_medium=referral&utm_content=GNiendorf/tracepy&utm_campaign=Badge_Grade_Dashboard)
[![Build Status](https://travis-ci.org/GNiendorf/tracepy.svg?branch=master)](https://travis-ci.org/GNiendorf/tracepy)
[![Documentation Status](https://readthedocs.org/projects/tracepy/badge/?version=latest)](https://tracepy.readthedocs.io/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/github/GNiendorf/tracepy/badge.svg?branch=master)](https://coveralls.io/github/GNiendorf/tracepy?branch=master)
[![PyPI version](https://badge.fury.io/py/tracepy.svg)](https://badge.fury.io/py/tracepy)

Ray Tracing and Optical Design in Python

## Overview

TracePy is a sequential ray tracing package written in Python 3 for designing optical systems in the geometric optics regime. It features lens optimization from Scipy.

## Installation

To use TracePy you can either clone the repository and use the command "pip3 install ." in the download directory, or you can download TracePy directly through pypi with the command below.

```
pip3 install tracepy
```

## Examples

To get started using the software, you can look at the examples provided in the 'examples' folder. The UI for TracePy is most similar to BEAM4, and TracePy's ray tracing algorithm was recreated from Spencer and Murty's iconic paper, "General Ray-Tracing Procedure".

## Contributing

I suggest reading Spencer and Murty's paper, ["General Ray-Tracing Procedure"](https://www.osapublishing.org/josa/viewmedia.cfm?uri=josa-52-6-672&seq=0) for an overview of TracePy's algorithm. There are also several open issues that we welcome collaborators on.
