# Frozen response source for the ER/NR band comparison

This is the **current local response snapshot used for the band study**, with
the corrected drift-time rejection condition. It is not an unmodified official
NEST or DIAMX release, and not LZ's production analysis. The exact copied source
hashes are in `source_copy_manifest.json`.

The wrapper's plotting bounds, corrected drift-time sampling, ER width-scale
switch, diagnostic ablation switch, and unused cascade interface are already
present in the frozen local source. This package does not introduce or activate
new physics changes. The public runner only calls `LZModel`, with scale 1 and
ablation 0. It does not call the cascade interface.

## Offline build and invocation

From the package root, on Linux or macOS with CMake >=3.16 and a C++17 compiler:

```
cmake -S response -B build_response
cmake --build build_response --parallel 2
./build_response/lz_response er 10000 response/inputs/flat_er_0p1_350.csv er_events.csv
./build_response/lz_response nr 10000 response/inputs/flat_nr_0p1_500.csv nr_events.csv
```

The source and all headers are bundled: no prebuilt library, WSL, Python
extension, network download, ROOT, Geant4, or installed NEST is required. WSL was
only the validation machine's Linux environment. Linux GNU 11.4.0/CMake 3.22.1
was actually tested; other compilers and platforms are not yet verified.
The unmodified bundled core can print many warnings for the NR input above
300 keVnr (the AmBe endpoint mentioned in its warning). Capture standard output
and standard error to a log when running millions of events. The NR 500 keVnr
reference was chosen to avoid truncation of the high-S1 conditional population;
it is not a new validation of NEST yields above that endpoint.

`lz_response er|nr EVENTS SPECTRUM.csv OUTPUT.csv` uses these fixed jobs:

| Population | NEST interaction | Seed | Included flat-energy input |
|---|---|---:|---|
| er | beta | 2452026 | 0.1--350 keV ER |
| nr | NR | 260903 | 0.1--500 keVnr |

Both call `LZModel("lz_ws2024", interaction, seed, spectrum, EVENTS, 1.0, 1.0, 0)`.
The selected detector calls `SetEFT2024Configuration()` in `LZ_WS2024.hh`.
Thus the historical filename is not a claim that the default detector settings
are used: inspect that function and the vectors it supplies to the wrapper.

## Output and filtering boundary

The output is one row per requested event, with 17 significant digits:

* `s1c_phd`: corrected S1 area in phd, exactly the wrapper's `scint[5]`.
* `s2c_phd`: corrected S2 in phd, exactly the wrapper's `scint2[7]`.
* `energy_rec`: the wrapper's original reconstructed energy in keV; this is
  **not true sampled energy**. Its S1 term uses `scint[7]` (spike count), so it
  must not be silently substituted for an area-based combined-energy estimate.
  It is not used to define or filter the fixed-S1 bands in this package.

The runner adds no absolute-value replacement, clipping, threshold, efficiency,
or selection. The copied wrapper's own conditions remain explicit in source:
`abs(scint[5]) > minS1 && scint[5] < maxS1` and the analogous S2 condition;
otherwise that observable is stored as -999. Negative subthreshold values can
also remain signed. The Python band estimator must select finite, positive
S1c and S2c before taking `log10(S2c)`; it should not turn negatives positive.
`analysis.hh` supplies minS1=1.5, minS2=42; the wrapper sets both upper bounds
to 1e7. No paper-specific event efficiency weights are supplied here.

The original ER precomputed export used six significant digits, while this
review runner preserves full double precision. A fresh full-precision band
estimate may consequently differ very slightly from the archived ER curve.
The root estimator's optional legacy rounding mode exists only to reproduce
that serialization boundary; it does not rescale physical fluctuations.

## Source and licensing boundary

`third_party/nest` contains the 9 exact core translation units and the complete
copied include tree from the locally used source, plus the original Apache-2.0
LICENSE, CITATION.md and upstream README. The local tree has no git metadata;
its file hashes, rather than an unverified tag label, define this snapshot.
`third_party/gcem` contains the locally used header-only dependency with its
original Apache-2.0 LICENSE and NOTICE. The wrapper is the local application
code; this README does not assign it an invented upstream release identity.

The response-only package contains no historical background spectra, no new1906
benchmark, and no O7/new-physics signal input. The flat files only define the
reference populations for the ER/NR band diagnostic.

## Fast review index

* `src/run_response.cpp`: the frozen job protocol and full precision CSV writer.
* `src/nestDiamxCore.cpp::LZModel`: energy sampling; positions and drift-time
  rejection; NEST GetYields/GetQuanta; detector GetS1/GetS2; signed output rules.
* `src/csvSpectra.hh::EnergySampler`: CSV loading, trapezoidal normalization,
  CDF construction and inverse-CDF interpolation (flat references here).
* `../third_party/nest/include/Detectors/LZ_WS2024.hh`:
  `SetEFT2024Configuration`, the 15/10/17-element NR/ER/width parameter vectors,
  detector gains/thresholds/maps, and uniform 96.5 V/cm drift field.
* `../third_party/nest/include/NEST/analysis.hh`: MCtruthPos=false,
  S1 Hybrid, S2 Full, original output thresholds.
* `../third_party/nest/src/NEST.cpp`: GetYields dispatches to GetYieldBetaGR or
  GetYieldNR; GetQuanta computes binomial + non-binomial recombination variance;
  GetS1 and GetS2 implement the detector response.
* `portable_validation.json`: independently compiled bundled source, ER/NR
  each 10,000 events, all three columns and complete CSV bytes identical to the
  previous clean-build reference (maximum numerical difference zero).
* `../data/small_reference/er_first100.csv` and `nr_first100.csv`: the first
  100 full precision rows from those verified jobs, for a fast seeded smoke
  comparison only. They are too small for the band result and are not used by
  the plot estimator. Their hashes are in `reference_manifest.json` beside them.
