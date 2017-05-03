SHMIP questionnaire
===================

Whole inter-comparison
----------------------

- *Name of the model & citation*
    Parallel Ice Sheet Model (PISM).

    .. code::

        @Article{Bueler.Pelt.2015,
          author    = {Bueler, E. and van Pelt, W.},
          title     = {Mass-conserving subglacial hydrology in the Parallel Ice
                       Sheet Model version 0.6},
          journal   = gmd,
          volume    = {8},
          number    = {6},
          pages     = {1613-1635},
          year      = {2015},
          doi       = {10.5194/gmd-8-1613-2015}
        }

- *Name of the modeller(s)*
    Julien Seguinot,
    Laboratory of Hydraulics, Hydrology and Glaciology,
    ETH Zürich, Switzerland.

- *On what computer was the model run*
    Vierzack04/05, ETH Zürich (Intel(R) Xeon(R) CPU E5-2687W v4 @ 3.00GHz)

- *Parallel or serial model*
    Parallel, on 8 cores.

- *Is the model available? If so, where?*
    Open-source, http://www.pism-docs.org/.

- *What revision/version of the model was used?*
    PISM stable 0.7.3.

- *Other remarks*
    In PISM, subglacial hydrology can be coupled to ice dynamics via a layer of
    porous, compressible till. The effective pressure felt by ice is the
    effective pressure on the till, which depends on water content in the till.
    At equilibrium, the till is typically saturated and the till effective
    pressure is a fixed fraction delta (default 0.02) of the overburden
    pressure. Instead, these results include effective pressure in the
    subglacial cavity system, although not that directly felt by ice.


Model tuning
------------

- *What model parameters were set as given in section Parameters?*
    All but the englacial void fracion of 1e-3 for sqrt and 1e-2 for valley.

- *Was the model tuned to A3 and/or A5? To N? q? Q?*
    No. PISM has similar physics to GlaDS except for the absence of channels,
    and yielded similar results for A3 without tuning.

- *What parameters where used for the tuning? Values?*
    None.

- *What are the values of any additional model parameters?*
    None.


Suite A
-------

- *How long was the model run for (model-time)?*
    5 years.

- *Are any parameters different from the ones used in Suite A?*
    Well, no.

- *CPU time used for the run.*
    - A1 0.69 h (x8 procs)
    - A2 1.68 h
    - A3 1.81 h
    - A4 0.80 h
    - A5 1.30 h
    - A6 8.71 h

- *How confident are you of model convergence?*
    Very confident for all runs.


Suite B
-------

- *How long was the model run for (model-time)?*
    5 years.

- *Are any parameters different from the ones used in Suite A?*
    No.

- *CPU time used for the run.*
    - B1 12.87 h (x8 procs)
    - B2 11.36 h
    - B3  9.05 h
    - B4  8.52 h
    - B5  2.15 h

- *How confident are you of model convergence?*
    Very confident for all runs.


Suite C
-------

- *How long was the model run for (model-time)?*
    30 days.

- *Are any parameters different from the ones used in Suite A?*
    No.

- *CPU time used for the run.*
    - C1 0.98 h (x8 procs)
    - C2 1.00 h
    - C3 0.99 h
    - C4 0.98 h

- *How confident are you of model convergence?*
    Very confident for C1 to C3; quite confident for C4.


Suite D
-------

- *How long was the model run for (model-time)?*
    5 years.

- *Are any parameters different from the ones used in Suite A?*
    No.

- *CPU time used for the run.*
    - D1 2.15 h (x8 procs)
    - D2 3.32 h
    - D3 4.69 h
    - D4 6.44 h
    - D5 8.38 h


- *How confident are you of model convergence?*
    Very confident for all runs.


Suite E
-------

- *How long was the model run for (model-time)?*
    5 years.

- *Are any parameters different from the ones used in Suite A?*
    Yes, the englacial void fracion of 1e-2.

- *Remarks*
    The high melt rate makes the model very slow. This is because it takes tiny
    time steps to respect the CFL condition on the last grid cell where
    effective diffusivity is very high, I assume. This is why I have used an
    increased englacial void fraction for valley runs.

- *CPU time used for the run.*
    - E1  12.24 h (x8 procs)
    - E2 109.89 h
    - E3  41.63 h
    - E4  52.98 h
    - E5  46.89 h


- *How confident are you of model convergence?*
    Very confident for all runs.


Suite D
-------

- *How long was the model run for (model-time)?*
    5 years.

- *Are any parameters different from the ones used in Suite A?*
    Yes, the englacial void fracion of 1e-2.

- *CPU time used for the run.*
    - F0 0.17 h (spin-up)
    - F1 0.98 h (x8 procs)
    - F2 2.07 h
    - F3 3.45 h
    - F4 5.02 h
    - F5 6.61 h

- *How confident are you of model convergence?*
    Very confident for all runs.
