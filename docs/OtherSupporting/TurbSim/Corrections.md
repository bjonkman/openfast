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

## User time-series input file: new optional `ExtrapFreqs` input

`ExtrapModel` fills in parts of the spectrum from a model that knows nothing about the measured
data, so the two generally do not join up. With `IECturbc = "10"` on a site whose measured
turbulence intensity was 1.7%, the IEC Kaimal spectrum was measured at **1100 times** the measured
spectrum for `u` and **3900 times** for `v` at the splice frequency.

`ExtrapFreqs` is an alternative for the components the file *does* supply: extrapolate the measured
spectra themselves rather than substituting a model. It takes two frequencies in Hz and may appear
either before or after `ExtrapModel`:

```
   0.10  0.30  ExtrapFreqs - Anchor and splice frequencies [Hz] for extrapolating the measured spectra
```

Above the **splice** frequency the measured spectra of components `1:nComp` are replaced by a power
law with the inertial-subrange slope (-5/3), whose level is set so that it reproduces the measured
band average around the **anchor** frequency. Setting the anchor to 0 (or omitting the line) leaves
the measured spectra alone, which is the previous behaviour.

`ExtrapModel` still handles the components the file does not supply, because there is no measured
spectrum to anchor those to. The two inputs are complementary and are normally used together.

### Choosing the two frequencies

They must come from the measured spectrum, not from preference:

* **anchor** — inside the band where the measurement actually follows -5/3.
* **splice** — the highest frequency you still trust.

For the nacelle-lidar data this was developed against, the measured slope was -1.66 and -1.65 over
0.05-0.15 Hz (physical), steepened to -2.4 and beyond above 0.15 Hz (instrument volume-averaging),
and reached -3.6 to -4.3 above 0.4 Hz (the roll-off of the interpolant used to remove the
instrument's zero-order hold). An anchor of 0.1 Hz and a splice of 0.3 Hz gave a step of 1.15 and
1.11 in `u` and `v`; leaving the splice at the 10 Hz Nyquist frequency of the file gave 45 and 122,
because at 10 Hz the "measured" spectrum being matched is numerical artefact rather than flow.

The slope is **imposed, not fitted**. Fitting it near the splice would extrapolate whichever
artefact dominates there. The trap that makes fitting unsafe in general: a zero-order-hold staircase
rolls off at roughly -1.8, close enough to -5/3 to look like an inertial subrange while being
meaningless, so a fitted slope would depend on how the input data happened to be resampled.

### Effect on exact reproduction of the input

A splice below the file's Nyquist frequency means the measured spectrum above it is discarded, so
**the input time series is no longer reproduced exactly** at the user-specified points. The phase
angles from the file are likewise used only up to the splice. This is usually a good trade, because
the discarded band is mostly artefact: in the case above, splicing at 0.3 Hz retained 0.1879 of the
0.1900 m/s measured standard deviation, i.e. it discarded 0.3% of the resolved energy. It is
nevertheless a change to what `TIMESR` normally guarantees, so it only happens when asked for.

The summary file records the splice and anchor frequencies and notes when exact reproduction has
been given up.

## Summary file: step in the spectrum at the splice frequency

The summary file now reports, per component at the reference point, the measured and extrapolated
spectral levels either side of the splice frequency and their ratio, and TurbSim issues a warning
when the ratio exceeds a factor of two either way:

```
  Spectrum step at the splice frequency (   0.300 Hz), reference point:
  (ratio is 1.0 if the spectrum continues smoothly with the -5/3 slope)
     component   measured        extrapolated       ratio
        u      7.8109E-03      7.4929E-03          1.15
        v      8.2502E-03      7.5830E-03          1.11
        w      2.0432E-01      1.7570E-01          1.04
```

Band averages are compared rather than single frequencies, because individual spectral estimates are
exponentially distributed and routinely differ from their mean by a factor of several. The ratio is
normalised by the step a smooth -5/3 spectrum would itself produce across the two averaging windows,
so 1.0 means continuous regardless of the splice frequency or window width. A component filled in
entirely by `ExtrapModel` therefore reports exactly 1.00, since both sides come from the same model.

A large ratio is not necessarily an error, but it always means one of two things: the level of the
fill-in spectrum does not match the site, or the splice frequency is above the highest frequency the
measurements genuinely resolve.
