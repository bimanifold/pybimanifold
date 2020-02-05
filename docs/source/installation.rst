Installation
============

Installation via Anaconda:
^^^^^^^^^^^^^^^^^^^^^^^^^^

Anaconda is a requirement for multiplex. You can install it `here <https://www.anaconda.com/distribution/>`_.

1) Create an empty Anaconda enviroment

.. code-block:: bash

   conda create --name myenv python=3.6.8 --no-default-packages

Python 3.6.8 is our testing version, but the code should work with any recent Python 3.x version.

2) Install the package inside the new enviroment:

.. code-block:: bash

  conda activate myenv
  conda install -c groever -c conda-forge multiplex

3) Test if the installation was successful

.. code-block:: python

   from multiplex import manifold as Manifold

Installation Problems? Consider Docker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`FEniCS <https://fenicsproject.org/>`_, the underlying CFD solver of multiplex, can be tricky to install depending on what OS and Python version you use. We test the successful installation of multiplex regularly on an Ubuntu 18.04 (Linux 64) operating
system with Python 3.6.8. You can recreate our testing environment with this DockerMake file.
