[![Build Status](https://travis-ci.org/bimanifold/pybimanifold.svg?branch=master)](https://travis-ci.org/bimanifold/pybimanifold)
[![codecov.io](https://codecov.io/gh/bimanifold/pybimanifold/branch/master/graph/badge.svg)](https://codecov.io/gh/bimanifold/pybimanifold)

# Multiplex

Multiplex is a CFD simulation package

## Requirements
OS: Linux 64 (Ubuntu 18.04 recommended)

## User information

0) Install [Anaconda](https://www.anaconda.com/distribution/ "Anaconda")

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

4) To test if the installation is successful, do:
```python
from multiplex import manifold as Manifold
```

## License
The license of this project is GNU Lesser General Public License v3.0, which is the license of [FEniCS](https://fenicsproject.org/), the main dependency of this project. However, we also give permission to distribute this source code under MIT License and BSD 3-Clause License.

## Developer information

Run dependency can be installed with:
```bash
conda install -c conda-forge fenics=2018.1.0 mshr=2018.1.0 matplotlib pyyaml
```

Build dependency can be installed with:
```bash
conda install -c conda-forge fenics=2018.1.0 mshr=2018.1.0 matplotlib pyyaml pytest pytest-cov codecov
```
