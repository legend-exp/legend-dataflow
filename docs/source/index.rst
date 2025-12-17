legend-dataflow
===============

*legend-dataflow* is a Python package based on Snakemake
`<https://snakemake.readthedocs.io/en/stable/index.html>`_ for running the data
production of LEGEND.  It is designed to calibrate and optimise hundreds of
channels in parallel before bringing them all together to process the data. It
takes as an input the metadata at `legend metadata
<https://github.com/legend-exp/legend-metadata>`_.

Getting started
---------------

Clone the repository using git.

.. code-block:: bash

   git clone https://github.com/legend-exp/legend-dataflow.git
   cd legend-dataflow

Then create a virtual environment to install *legend-dataflow* to.
Use e.g. ``uv`` for that:

.. code-block:: bash

   uv venv --python 3.12
   source .venv/bin/activate
   uv pip install -e ".[dev]"

With ``[dev]`` you install the development dependencies. You might want to
use a different set of dependencies dependent on you use case.

Adapt the ``dataflow-config.yaml`` and add a workflow profile in
``workflow/profiles/`` if you want to set the dataflow up for a new host.
Otherwise, check if your host is already configured or if ``bare`` applies for you.

Install the dataflow using

.. code-block:: bash

   dataflow -v install -s <host> dataflow-config.yaml

with ``<host>`` being the hostname as configured in ``dataflow-config.yaml``.
This command installs all the necessary software to run the dataflow to
``.snakemake/legend-dataflow/venv``.
Be sure to clear the numba cache (defined in the config) in case of software updates.



Next steps
----------

.. toctree::
   :maxdepth: 1

   Package API reference <api/modules>

.. toctree::
   :maxdepth: 1

   user_manual

.. toctree::
   :maxdepth: 1
   :caption: Related projects

   LEGEND Data Objects <https://legend-pydataobj.readthedocs.io>
   Decoding Digitizer Data <https://legend-daq2lh5.readthedocs.io>
   Digital Signal Processing <https://dspeed.readthedocs.io>
   Pygama <https://pygama.readthedocs.io>

.. toctree::
   :maxdepth: 1
   :caption: Development

   developer
   Source Code <https://github.com/legend-exp/legend-dataflow>
