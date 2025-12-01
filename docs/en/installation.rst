Installation
============

Requirements
------------

* Python >= 3.10
* ORCA 6.0.1 installed and available in PATH
* RDKit >= 2023.0.0

Installing ORCA
----------------

ORCA is a quantum chemistry program that must be installed separately. 

1. Download ORCA from the official website: https://orcaforum.kofo.mpg.de/
2. Extract and add the ORCA executable to your system PATH
3. Verify installation by running::

   orca --version

Installing the Library
-----------------------

Using pip
~~~~~~~~~

.. code-block:: bash

   pip install orca-descriptors

Using Poetry (Development)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone <repository-url>
   cd orca_descriptors
   poetry install

Verification
------------

After installation, verify that the library works::

   python -c "from orca_descriptors import Orca; print('OK')"

You can also test the CLI::

   orca_descriptors --help

