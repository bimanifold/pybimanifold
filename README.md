[![Build Status](https://travis-ci.org/bimanifold/pybimanifold.svg?branch=master)](https://travis-ci.org/bimanifold/pybimanifold)
[![codecov.io](https://codecov.io/gh/bimanifold/pybimanifold/branch/master/graph/badge.svg)](https://codecov.io/gh/bimanifold/pybimanifold)

This is a public repo for multiplexing microfluidics droplet makers based on flow flocusing. Please reach out with any questions.

author: Benedikt Groever

email: -

1) Create an empty Anaconda enviroment

```bash
conda create --name multiplex python=3.6.8 --no-default-packages
```

2) Activate Anaconda enviroment

```bash
conda activate multiplex
```

3) You can install the package as:
```bash
conda install -c groever -c conda-forge multiplex
```

4) Test if the build is successful with:
```python
from multiplex import manifold as Manifold
```

Additional info:
Run dependency can be installed with
```bash
conda install -c conda-forge fenics=2018.1.0 mshr=2018.1.0 matplotlib pyyaml
```

Build dependency can be installed with
```bash
conda install -c conda-forge fenics=2018.1.0 mshr=2018.1.0 matplotlib pyyaml pytest pytest-cov codecov
```
