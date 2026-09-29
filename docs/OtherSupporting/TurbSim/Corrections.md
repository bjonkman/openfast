# Corrections to TurbSim_v2.00.pdf DRAFT

## p. 34: Input file for User-Defined Profiles
The User-Defined profiles listed in the "Profiles" section contain horizontal angles, not wind direction, so it should say
  `3. Horizontal wind angle (degrees, measured counter-clockwise from above)`

## User time-series input file: new optional `ExtrapModel` input

The user time-series input file (`TurbModel = "TIMESR"`) takes two optional inputs, which may appear
in either order on the lines after `RefPtID`:

```
   "IECKAI"  ExtrapModel - Spectral model used for the velocity components this file does not supply ("NONE","IECKAI","IECVKM","SMOOTH")
0.1000 0.3000  ExtrapFreqs - Anchor and splice frequencies [Hz] for extrapolating the measured spectra of the components this file does supply
```

They divide the work as follows, and are normally used together:

| Part of the spectrum | Filled by |
|---|---|
| components above `nComp` (`w` when `nComp` = 2; `v` and `w` when `nComp` = 1), all frequencies | `ExtrapModel` |
| components `1:nComp` above the splice frequency, when `ExtrapFreqs` is given | the -5/3 fit (see below) |
| components `1:nComp` above the file's Nyquist frequency, when `ExtrapFreqs` is **not** given | `ExtrapModel` |

Both are optional. Omit them (or set `ExtrapModel = "NONE"` and the anchor frequency to 0) and the
unspecified parts of the spectrum contain no energy, so files written for the original format are read
and produce the same results as before.

`IECKAI` and `IECVKM` take their turbulence level from the `IECstandard` and `IECturbc` inputs of the
primary input file, exactly as in a normal IEC run, so `IECturbc` sets how much turbulence is added.
(Those inputs are otherwise ignored for `TIMESR`.) `SMOOTH` instead uses `RICH_NO` and `UStar`. The
coherence of the run is unaffected: it continues to come from the `SCMod`/`InCDec` inputs for the
`TIMESR` model rather than from the fill-in model.

`API` is **not** accepted. `Spec_API` is written against a 10 m reference height and a 1-hour mean
`URef`, the checks that enforce those conditions test `p%met%TurbModel_ID` and so never fire when API
is only the fill-in model, and `TIMESR` has by then overwritten `RefHt` and `URef` with the reference
point's height and mean speed. The resulting spectrum would be scaled by an arbitrary factor.

### Reaching the frequencies above the file's Nyquist

Run with `TimeStep` **smaller** than the time between rows of the time-series file. The file's time
step must be an integer multiple of `TimeStep`, and `AnalysisTime` must give both series the same
duration so that they share a frequency increment. For example, a file written with a 0.05 s step and
a run with `TimeStep = 0.0125` resolves the field to 40 Hz, so the 10-40 Hz band gets filled in.
Previously `TimeStep` was required to equal the file's time step, which made this band unreachable.

Reducing `TimeStep` is **not** required in order to use `ExtrapFreqs`. A splice frequency below the
file's own Nyquist frequency already creates a band for the -5/3 fit to occupy, so `ExtrapFreqs` is
useful at `TimeStep` equal to the file's time step.

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
law with the inertial-subrange slope (-5/3). Its level is fitted so that it reproduces the measured
arithmetic mean over the band **from the anchor frequency up to the splice frequency**. Setting the
anchor to 0 (or omitting the line) leaves the measured spectra alone, which is the previous
behaviour.

The fit band deliberately *ends* at the splice. A band that stops short of it leaves the power law
free to drift away from the measurements over the intervening gap, and because these spectra roll off
more steeply than -5/3 the drift is always upward. Measured with anchor 0.1 Hz and splice 0.3 Hz,
fitting over a band centred on the anchor (0.067 to 0.15 Hz) left the extrapolation 1.32x high at the
splice; fitting over 0.1 to 0.3 Hz brought that to 1.14x.

`ExtrapModel` still handles the components the file does not supply, because there is no measured
spectrum to anchor those to. The two inputs are complementary and are normally used together.

### Choosing the two frequencies

Together they define the fit band, so they must come from the measured spectrum, not from preference:

* **anchor** — the lowest frequency you consider clean and representative.
* **splice** — the highest frequency you still trust. Also the top of the fit band.

For the nacelle-lidar data this was developed against, the spectrum is closest to -5/3 somewhere
around 0.05 to 0.15 Hz, steepens past -2 above that (instrument volume-averaging), and steepens
further above roughly 0.4 Hz (the roll-off of whatever interpolant removed the instrument's
zero-order hold). Note the measured slope depends on the band edges you pick *and* on how the hold
was removed: over 0.05-0.15 Hz it was -1.65, but over 0.067-0.15 Hz already -1.84, and with the raw
staircase retained instead of interpolated, -1.92. Establish it for your own data rather than
assuming these numbers.

An anchor of 0.1 Hz with a splice of 0.3 Hz is a reasonable starting point for this instrument.
Leaving the splice at the file's Nyquist frequency is not: the "measured" spectrum being matched
there is numerical artefact rather than flow, and the step reached factors of 45 and 122 in `u` and
`v`.

Expect the power law to remain somewhat above the measured continuation well above the splice — about
2.6x by 5 Hz in the case above. That is a *slope* mismatch, not a level error, so it cannot be fixed
by moving the fit band, and it is not necessarily wrong: the measured continuation up there is the
instrument's attenuation plus the output-format quantisation floor, which is what the extrapolation
exists to replace.

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
