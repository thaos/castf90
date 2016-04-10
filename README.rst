=======
castf90
=======

Circulation analogue simulation tool in fortran95. 

|Install with Conda|

Installing castf90
==================

Anaconda
--------

castf90 is available as conda package. Install it with the following command:

.. code:: bash

    $ conda install -c birdhouse -c ioos castf90

From github
-----------

Clone the castf90 github repo:

.. code:: bash

    $ git clone https://github.com/sradanov/castf90.git
    $ cd castf90

Prepare a conda environment with the dependencies and activate it:

.. code:: bash

    $ conda env create -f environment.yml
    $ source activate castf90

or install the castf90 dependencies manually in your current conda enviroment:

.. code:: bash

   $ conda install -c birdhouse gcc openblas netcdf-fortran lapack95 cdo nco
   $ conda install -c ioos cdo nco

Build the analogue.out executable:

.. code:: bash

   $ make -f Makefile.conda


.. |Install with Conda| image:: https://anaconda.org/birdhouse/castf90/badges/installer/conda.svg
   :target: https://anaconda.org/birdhouse/castf90
