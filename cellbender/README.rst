.. _installation:

CellBender Installation
============

Manual Installation
-------------------

The recommended installation is as follows. Create a conda environment and activate it:

.. code-block:: console

  $ conda create -n CellBender python=3.7
  $ source activate CellBender

Install the `pytables <https://www.pytables.org>`_ module:

.. code-block:: console

  (CellBender) $ conda install -c anaconda pytables

Install `pytorch <https://pytorch.org>`_ (shown below for CPU; if you have a CUDA-ready GPU, please skip
this part and follow `these instructions <https://pytorch.org/get-started/locally/>`_ instead):

.. code-block:: console

   (CellBender) $ conda install pytorch torchvision -c pytorch

Clone this repository and install CellBender:

.. code-block:: console

   (CellBender) $ git clone https://github.com/broadinstitute/CellBender.git
   (CellBender) $ pip install -e CellBender

Running CellBender ``remove-background`` on Placenta Dataset
--------------

First, download 10x Genomics human placenta datasets:

.. code-block:: console

   (CellBender) $ cd ../data/; bash download_placenta.sh

Then, you can run CellBender on all six replicates with the following command:

.. code-block:: console

   (CellBender) $ cd ../cellbender; bash cellbender_run.sh

**NOTE:** Each CellBender command takes a while to run. Plan on these analyses taking up to 12 hours for all 6 replicates.

Converting CellBender outputs to AnnData for downstream analysis
--------------

Finally, you'll need to put the ``.h5`` outputs into ``.h5ad`` format in the ``data/`` directory for ease of ``dropkick`` filtering:

.. code-block:: console

   (CellBender) $ python h5_to_h5ad.py; mv *.h5ad ../data/
