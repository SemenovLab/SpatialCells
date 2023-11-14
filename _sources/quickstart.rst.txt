Quick Start Guide
=================

Installation
------------

SpatialCells can currently only be installed from source. 
To install, download the repository and run the following commands from the root directory of SpatialCells:

.. code-block::

    pip install -r requirements.txt
    pip install .


It is recommended to install SpatialCells in a virtual environment, such as `conda <https://docs.conda.io/en/latest/>`_.
A conda environment can be created with the following yaml file:

.. code-block::

    conda env create -f conda.yml

The conda.yml is specified as follows:

.. literalinclude:: _static/conda.yaml
  :language: yaml

