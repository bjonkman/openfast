# Corrections to TurbSim_v2.00.pdf DRAFT

## p. 34: Input file for User-Defined Profiles
The User-Defined profiles listed in the "Profiles" section contain horizontal angles, not wind direction, so it should say
  `3. Horizontal wind angle (degrees, measured counter-clockwise from above)`

## User time-series input file: new optional `ExtrapModel` input

The user time-series input file (`TurbModel = "TIMESR"`) takes an optional input on the line
after `RefPtID`:

```
     "NONE"  ExtrapModel - Spectral model used to fill in the velocity components and frequencies not specified in this file ("NONE","IECKAI","IECVKM","API","SMOOTH")
```

It names the spectral model used to fill in the parts of the spectrum the file does not specify:

1. velocity components above `nComp` (i.e. `w` when `nComp` = 2, and `v` and `w` when `nComp` = 1),
   at all frequencies; and
2. all components at frequencies above the Nyquist frequency of the time series in the file.

The input is optional. When the line is omitted the model defaults to `NONE` and the unspecified
parts of the spectrum contain no energy, so files written for the original format are read and
produce the same results as before.

To reach case 2, run with `TimeStep` **smaller** than the time between rows of the time-series
file. The file's time step must be an integer multiple of `TimeStep`, and `AnalysisTime` must give
both series the same duration so that they share a frequency increment. For example, a file
written with a 0.05 s step and a run with `TimeStep = 0.0125` resolves the field to 40 Hz and
synthesises the 10-40 Hz band. Previously `TimeStep` was required to equal the file's time step,
which made this case unreachable.

`IECKAI`, `IECVKM` and `API` take their turbulence level from the `IECstandard` and `IECturbc`
inputs of the primary input file, exactly as in a normal IEC run, so `IECturbc` sets how much
turbulence is added. (Those inputs are otherwise ignored for `TIMESR`.) `SMOOTH` instead uses
`RICH_NO` and `UStar`. The coherence of the run is unaffected: it continues to come from the
`SCMod`/`InCDec` inputs for the `TIMESR` model rather than from the fill-in model.

Note that this feature was present in the source but inert before this change: the model was
read from `p%usr%TurbModel_ID`, which was initialised to `SpecModel_NONE` and never assigned, so
`Spec_TimeSer_Extrap` always returned a zero spectrum. This is why the `nComp` comment in the
file header ("if < 3 other components will be generated using values from input file") did not
previously describe the observed behavior.
