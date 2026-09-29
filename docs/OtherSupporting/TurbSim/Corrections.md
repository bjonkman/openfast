# Corrections to TurbSim_v2.00.pdf DRAFT

## p. 34: Input file for User-Defined Profiles
The User-Defined profiles listed in the "Profiles" section contain horizontal angles, not wind direction, so it should say
  `3. Horizontal wind angle (degrees, measured counter-clockwise from above)`

## User time-series input file: new optional `ExtrapModel` input

The user time-series input file (`TurbModel = "TIMESR"`) takes two optional inputs, which may appear
on the lines after `RefPtID`:

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
`TIMESR` model rather than from the fill-in model. `API` is **not** accepted. 

### Reaching the frequencies above the file's Nyquist

Run with `TimeStep` **smaller** than the time between rows of the time-series file. The file's time
step must be an integer multiple of `TimeStep`, and `AnalysisTime` must give both series the same
duration so that they share a frequency increment. For example, a file written with a 0.05 s step and
a run with `TimeStep = 0.0125` resolves the field to 40 Hz, so the 10-40 Hz band gets filled in.
Previously `TimeStep` was required to equal the file's time step, which made this band unreachable.

Reducing `TimeStep` is **not** required in order to use `ExtrapFreqs`. A splice frequency below the
file's own Nyquist frequency already creates a band for the -5/3 fit to occupy, so `ExtrapFreqs` is
useful at `TimeStep` equal to the file's time step.

## User time-series input file: new optional `ExtrapFreqs` input

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
more steeply than -5/3 the drift is always upward. 

`ExtrapModel` still handles the components the file does not supply, because there is no measured
spectrum to anchor those to. The two inputs are complementary and are normally used together.

### Choosing the two frequencies

Together they define the fit band, so they must come from the measured spectrum, not from preference:

* **anchor** — the lowest frequency you consider clean and representative.
* **splice** — the highest frequency you still trust. Also the top of the fit band.

Note that measured spectra may be affected by instrument volume-averaging and processing (e.g., removing zero-order hold).


### Effect on exact reproduction of the input

A splice below the file's Nyquist frequency means the measured spectrum above it is discarded, so
**the input time series is no longer reproduced exactly** at the user-specified points. The phase
angles from the file are likewise used only up to the splice. This is usually a good trade, because
the discarded band is mostly artefact. It is a change to what `TIMESR` normally guarantees, so it 
only happens when asked for.

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

A large ratio is not necessarily an error, but it always means one of two things: the level of the
fill-in spectrum does not match the site, or the splice frequency is above the highest frequency the
measurements genuinely resolve.

It is recommended that users calculate and plot the PSDs of the resulting time series when using this
feature.
