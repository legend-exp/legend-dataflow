Pipeline Overview
=================

.. contents:: Contents
   :local:
   :depth: 2


Introduction
------------

The legend-dataflow pipeline transforms raw digitiser data into physics-ready event
data through a sequence of processing stages called *tiers*. Each tier builds on the
previous one, progressively applying calibrations and higher-level reconstruction.

The pipeline has two interleaved tracks:

1. **Calibration parameter generation** – derives per-channel calibration parameters
   from dedicated calibration runs (``cal`` datatype).
2. **Data processing** – applies those parameters to both calibration and physics data
   (``phy`` datatype).

For each tier, parameter generation must complete before the corresponding data
processing step. Parameters are stored in YAML files alongside validity windows that
record which runs they apply to.


Data Format
-----------

All tier data are stored in `LH5 <https://legend-pydataobj.readthedocs.io>`_ (LEGEND
HDF5) format. Within each file, data are organised by channel using the table naming
convention specified in ``table_format`` in the configuration, for example:

- ``ch{ch:07d}/raw`` – raw waveform data for channel ``ch``
- ``ch{ch:07d}/dsp`` – DSP-processed data for channel ``ch``
- ``{grp}/evt`` – event-level data for detector group ``grp``

File names follow the pattern::

   {experiment}-{period}-{run}-{datatype}-{timestamp}-tier_{tier}.lh5


File Key Structure
------------------

Files are identified by *file keys* composed of:

- **experiment** – the experiment name (e.g. ``l200``)
- **period** – data-taking period (e.g. ``p03``)
- **run** – run number within the period (e.g. ``r001``)
- **datatype** – ``cal`` (calibration) or ``phy`` (physics)
- **timestamp** – UTC timestamp of the DAQ cycle (``YYYYMMDDTHHMMSSz`` format)

The Snakemake target label (e.g. ``all-l200-p03-r001-phy``) selects a set of file
keys via the runlist defined in the ``legend-datasets`` repository.


Calibration Groupings and Partitions
--------------------------------------

Some calibration parameters are more stable when derived from multiple runs together
rather than a single run. The *partition* concept groups runs that share the same
calibration parameters (e.g. because detector conditions were stable across those
runs).

Partition assignments are defined in ``cal_groupings.yaml`` in the ``legend-datasets``
repository and loaded by the ``CalGrouping`` class. Rules with ``part`` in their name
(e.g. ``psp``, ``pht``) operate at partition level.


Processing Tiers
----------------

RAW tier
~~~~~~~~~

**Rule file:** ``rules/raw.smk``

Converts raw DAQ files to LH5 format:

- **ORCA format** (``build_raw_orca``) – decodes data from the ORCA DAQ system used
  for HPGe detectors.
- **FCIO format** (``build_raw_fcio``) – decodes data from the FCIO DAQ system used
  for SiPM detectors.
- **Blinding** (``build_raw_blind``) – applies an energy-dependent random shift to
  high-energy events in physics data to blind the analysis region of interest. The
  blinding calibration coefficients are derived separately (see
  ``rules/blinding_calibration.smk``).

Each raw file preserves the waveform data and DAQ metadata for all channels recorded
in that DAQ cycle.

TCM tier
~~~~~~~~~

**Rule file:** ``rules/tcm.smk``

Builds the *Timing Correction Module* (TCM), which provides cross-channel timing
alignment information. It also identifies pulser events (artificial periodic signals
injected for monitoring) via ``build_pulser_ids``, which are used to track detector
stability throughout processing.

DSP tier
~~~~~~~~~

**Rule files:** ``rules/dsp.smk``, ``rules/dsp_pars_geds.smk``,
``rules/dsp_pars_spms.smk``

*Parameter generation* (per-channel, from calibration data):

- **HPGe detectors** – derives optimal digital filter parameters for energy
  reconstruction: pole-zero correction, energy filter shaping, charge trapping
  correction. Each channel is optimised independently.
- **SiPM detectors** – derives trigger threshold parameters optimised per channel.

*Data processing:*

Applies the per-channel DSP parameters to all raw data, producing DSP-tier files
containing quantities such as:

- Calibrated energy estimates
- Baseline and pile-up flags
- Waveform shape parameters

PSP tier
~~~~~~~~~

**Rule files:** ``rules/psp.smk``, ``rules/psp_pars_geds.smk``

The *Partition-level DSP* (PSP) tier re-derives DSP parameters by averaging over all
calibration runs within a partition. This improves parameter stability for detectors
with sufficient statistics across multiple runs. The averaged parameters are then
applied to data to produce PSP-tier files with the same structure as DSP.

HIT tier
~~~~~~~~~

**Rule files:** ``rules/hit.smk``, ``rules/hit_pars_geds.smk``

*Parameter generation* (per-channel, from calibration data):

- **Energy calibration** – fits gamma-ray peaks to establish the conversion from ADC
  counts to energy in keV.
- **Pulse shape discrimination (PSD)** – derives parameters for separating signal-like
  events from backgrounds based on waveform shape.

*Data processing:*

Applies calibration and PSD parameters, producing HIT-tier files with calibrated
energies and PSD quantities for each detector channel.

PHT tier
~~~~~~~~~

**Rule files:** ``rules/pht.smk``, ``rules/pht_pars_geds.smk``,
``rules/pht_pars_geds_fast.smk``

The *Partition-level HIT* (PHT) tier is the most complex calibration stage. It derives
advanced analysis parameters at partition level:

- **Quality control** (``qc``) – channel-level data quality flags
- **Physics QC** (``qc_phy``) – event-level quality cuts for physics analysis
- **Partition energy calibration** (``ecal_part``) – refined energy calibration using
  the full partition statistics
- **Alpha-over-energy** (``aoe``) – a pulse shape parameter used to discriminate
  alpha-particle surface events
- **Log-quality** (``lq``) – a complementary pulse shape parameter

A fast mode (``pht_pars_geds_fast.smk``) is available for quicker turnaround at
reduced accuracy.

ANN tier
~~~~~~~~~

**Rule file:** ``rules/ann.smk``

Applies artificial neural network (ANN) cuts to data from coaxial HPGe detectors.
These cuts provide additional background rejection beyond classical PSD. The ANN
produces both channel-level (``ann``) and partition-level (``pan``) outputs.

EVT tier
~~~~~~~~~

**Rule file:** ``rules/evt.smk``

Builds *event-level* data by combining information across all detector channels:

- **Cross-talk correction** – removes electrical cross-talk between channels
- **Multi-detector coincidences** – identifies events where energy is deposited in
  multiple detectors simultaneously
- **Event energy** – computes the total energy deposited in each event

The EVT tier produces both run-level (``evt``) and partition-level (``pet``) outputs.

SKM tier
~~~~~~~~~

**Rule file:** ``rules/skm.smk``

The *skim* tier applies the final physics event selection, keeping only events that
pass all quality and analysis cuts. This produces compact files suitable for
physics analysis.


Parameter Catalogs
------------------

Calibration parameters are tracked through *parameter catalogs* – YAML files that
record which parameter set is valid for which time range. At the start of each
Snakemake run, validity catalog files are written to:

- ``{par_path}/dsp/validity.yaml``
- ``{par_path}/hit/validity.yaml``
- ``{par_path}/psp/validity.yaml``
- ``{par_path}/pht/validity.yaml``

These catalogs are used during data processing to look up the correct parameters for
each file's timestamp. The ``ParsKeyResolve`` class (in
``methods/create_pars_keylist.py``) and ``ParsCatalog`` class (in
``methods/pars_loading.py``) implement catalog resolution.


Channel Merging
---------------

**Rule file:** ``rules/channel_merge.smk``

After per-channel parameter generation, individual channel parameter files and
diagnostic plots are merged into single combined files per tier. This simplifies
downstream handling and allows the full detector array's parameters to be reviewed
in one place.


Blinding
--------

**Rule files:** ``rules/blinding_calibration.smk``, ``rules/blinding_check.smk``

To prevent analysis bias, events in the ``0νββ`` signal region of interest are
energy-blinded before physics analysis. The blinding procedure:

1. Derives a calibration curve from ``cal`` data (``blinding_calibration.smk``)
2. Validates the curve (``blinding_check.smk``)
3. Applies a reproducible energy-dependent shift to candidate signal events in
   ``phy`` data during raw-tier conversion (``raw.smk``)

Blinding is applied early in the pipeline so that all downstream tiers see only
blinded data.


Processing Flow Summary
-----------------------

The following diagram summarises the dependency structure of the main processing stages:

.. code-block:: text

                    DAQ files (ORCA / FCIO)
                            |
                     [raw.smk] RAW
                            |
                     [tcm.smk] TCM
                       /         \
      [dsp_pars_geds/spms]    DSP params
                       \         /
                     [dsp.smk] DSP
                       /         \
            [psp_pars_geds]    PSP params
                       \         /
                     [psp.smk] PSP
                       /         \
            [hit_pars_geds]    HIT params
                       \         /
                     [hit.smk] HIT
                       /         \
            [pht_pars_geds]    PHT params
                       \         /
                     [pht.smk] PHT
                            |
                     [ann.smk] ANN
                            |
                     [evt.smk] EVT
                            |
                     [skm.smk] SKM
