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

    $ conda install -c birdhouse -c conda-forge castf90

From github
-----------

Clone the castf90 github repo:

.. code:: bash

    #$ git clone https://github.com/sradanov/castf90.git
    $ git clone https://github.com/bird-house/castf90.git
    $ git checkout pingudev
    $ cd castf90

Create the conda environment `castf90` with necessary dependencies and activate it:

.. code:: bash

    $ conda env create -f environment.yml
    $ source activate castf90

or install the `castf90` dependencies manually in your current conda enviroment:

.. code:: bash

   $ conda install -c birdhouse -conda-forge gcc mpl openblas netcdf-fortran lapack95 cdo nco

Build the analogue.out executable:

.. code:: bash

   $ make -f Makefile.conda clean
   $ make -f Makefile.conda          # using intel MPL
   or
   $ make -f Makefile.conda NOMKL=1  # using openblas


.. |Install with Conda| image:: https://anaconda.org/birdhouse/castf90/badges/installer/conda.svg
   :target: https://anaconda.org/birdhouse/castf90
