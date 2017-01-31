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
    Piz Daint, Swiss National Supercomputing Centre.

- *Parallel or serial model*
    Parallel, on 36 core (probably too many).

- *Is the model available? If so, where?*
    Open-source, http://www.pism-docs.org/.

- *What revision/version of the model was used?*
    PISM stable 0.7.3.

- *Other remarks*
    In PISM, subglacial hydrology can be coupled to ice dynamics via a layer of
    porous, compressible till. These experiments were conducted without a till
    layer, but in the future I would like to submit a second set of experiments
    including a till layer.


Model tuning
------------

- *What model parameters were set as given in section Parameters?*
    All but the englacial void fracion which was set to 1e-3.

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
    - A1 0.39 h (x36 procs)
    - A2 0.97 h
    - A3 1.04 h
    - A4 0.45 h
    - A5 0.76 h
    - A6 5.21 h

- *How confident are you of model convergence?*
    Very confident for all runs.


Suite B
-------

- *How long was the model run for (model-time)?*
    5 years.

- *Are any parameters different from the ones used in Suite A?*
    No.

- *CPU time used for the run.*
    - B1 7.38 h (x36 procs)
    - B2 6.53 h
    - B3 5.14 h
    - B4 4.95 h
    - B5 1.19 h

- *How confident are you of model convergence?*
    Very confident for all runs.


Suite E
-------

- *How long was the model run for (model-time)?*
    - E1 1.3 a
    - E2 0.5 a
    - E3 0.5 a
    - E4 0.2 a
    - E5 0.3 a

- *Are any parameters different from the ones used in Suite A?*
    No.

- *Remarks*
    The high melt rate makes the model very slow. This is because it takes tiny
    time steps to respect the CFL condition on the last grid cell where and
    Effective diffusivity is very high, I assume. If I have a chance I will try
    to improve this in the future, perhaps increasing porosity.

- *CPU time used for the run.*
    - E1 23.55 h (x36 procs)
    - E2 21.03 h
    - E3 22.81 h
    - E4 20.13 h
    - E5 21.08 h

- *How confident are you of model convergence?*
    Quite confident for E1, E2 and E3.
    Not confident for E4 and E5.
