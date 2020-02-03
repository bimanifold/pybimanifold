Installation
============

Installation via Anaconda:
^^^^^^^^^^^^^^^^^^^^^^^^^^

Anaconda is a requirement for this package. You can install it `here <https://www.anaconda.com/distribution/>`_.

1) Create an empty Anaconda enviroment

.. code-block:: bash

   conda create --name multiplex python=3.6.8 --no-default-packages

2) Install the package inside the new enviroment:

.. code-block:: bash

  conda activate multiplex
  conda install -c groever -c conda-forge multiplex

3) Test if the installation was successful

.. code-block:: python

   from multiplex import manifold as Manifold

Installation via Docker
^^^^^^^^^^^^^^^^^^^^^^^

FEniCS, the underlying CFD solver is a large package, which can be complicated
to install. We test our package regularly on an Ubuntu 18.04 (Linux 64) host operating
system. You can recreate our testing environment using Docker.
