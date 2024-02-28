Quick Start Guide
=================

Installation
------------

SpatialCells can currently only be installed from source. 
To install, download code from the `repository <https://github.com/SemenovLab/SpatialCells>`_ and change to the code folder.

Run the following commands from the root directory of SpatialCells:

.. code-block::

    pip install -r requirements.txt
    pip install .


It is recommended to install SpatialCells in a virtual environment, such as `conda <https://docs.conda.io/en/latest/>`_.
A conda environment can be created with the following yaml file:

.. code-block::

    conda env create --name spatialcells --file=conda.yaml
    pip install .

The conda.yaml is specified as follows:

.. literalinclude:: ../conda.yaml
  :language: yaml

