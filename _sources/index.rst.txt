.. spatialcells documentation master file, created by
   sphinx-quickstart on Thu Nov  9 16:10:17 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SpatialCells
========================================

SpatialCells is an open-source software package designed to perform region-based exploratory 
analysis and characterization of tumor microenvironments using multiplexed single-cell data. The package uses 
the `Anndata <https://anndata.readthedocs.io/en/latest/>`_ framework and complements existing 
single-cell analysis packages such as `MCMICRO <https://mcmicro.org/>`_ and `SCIMAP <https://scimap.xyz/>`_. 
In particular, SpatialCells can efficiently analyze samples comprising millions of cells 
without manual annotation and automatically extract quantitative features of cell regions, 
enabling subsequent association analyses and machine learning predictions at scale.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   tutorials
   spatialcells



* :ref:`genindex`