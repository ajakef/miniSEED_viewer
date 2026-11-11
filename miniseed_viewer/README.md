# miniseed_viewer

A quick-look viewer for miniSEED files, in the spirit of PASSCAL's PQL.

    miniseed_viewer file1.mseed file2.mseed ...

Files can also be dropped onto the window. Opening files only reads their headers (fast); the samples of
the channels on the current page are read when the page is shown, and old pages are dropped again, so
large sets of files start quickly. "Traces per screen" (Settings) sets the page size; use the
Previous/Next buttons (or PageUp/PageDown) to move between pages.

## Using it

| Tab | What it does |
|---|---|
| Main | Raw traces. Click once to set t1, click again to set t2; the other tabs then show that window as soon as you open them. Unchecked traces are faded and left out of the other tabs. "Calibration..." sets a bitweight (physical units per count) and unit name per channel; every tab then shows physical units. |
| Analyses | Makes new results from two of the traces ticked on Main. A **Filter** (any saved filter, with "Edit filters...") is chosen first: it is applied to the trace plots, the cross-correlation and the correlogram (the whole record is filtered, then the interval is cut out, so start-up effects stay outside), and for the difference only to find R (the difference itself is never filtered). Choose trace 1 and trace 2 (order matters), optionally give traces a time shift (positive = later), optionally limit to the window selected on Main (in "seconds since start" mode the window is converted to clock time with trace 1's start and applied to both traces). **Difference** adds a new trace DIFFn = trace 1 - trace 2 (if the sample rates differ, the faster trace is low-pass filtered at 40% of the slower rate and resampled to it; trace 2 is interpolated onto trace 1's sample times); with "Find the scale factor R" ticked, the output is sqrt(R)·trace 1 − trace 2/sqrt(R), where R = std(trace 2)/std(trace 1) over the shared time is the value that minimizes the output's RMS after its mean is removed (a saved filter can be chosen to be applied to both traces for finding R only; the output is always the unfiltered scaled difference) (R is shown in the message and the trace's description); it appears on Main and can be shown in Magnify, Spectra, Spectrogram and Split. **Cross-correlation** plots c(lag) = sum of x1(t)·x2(t+lag) (means removed) for the lag range you give, fitted to the screen; a peak at a positive lag means trace 2 is a delayed copy of trace 1. By default only the time interval currently shown in the two trace plots is correlated (zoom there to choose it; untick the box to use the whole traces); the interval used is reported next to the result. The two chosen traces are plotted below the controls with their time shifts applied (calibrated, polarity-corrected); click twice on a plot to zoom to the interval between the clicks (as in Magnify), or use the mouse wheel; right-drag pans; zooming affects the time axis only, each trace's vertical scale always fits what is on screen, "Show all" and "Zoom to the window selected on Main" set the view. **Correlogram** is like a spectrogram of cross-correlations: time along the same axis as the trace plots, time lag up the side, and the correlation coefficient (-1 to 1, red-blue) as color. Parameters: window length (s), overlap (%), smallest and largest lag; it uses the time interval shown in the trace plots (zoom first, then Calculate). The window should be much longer than the largest lag; windows that are mostly empty are left blank. The color bar autoscales by default (symmetric about 0, to the 99.5th percentile of |coefficient|, so a few strong peaks saturate rather than flattening everything else); untick "autoscale color bar" to type a limit. Results are not saved between runs. |
| Magnify | Filtering, detrending and polarity reversal (for all traces, or per trace with each trace's "settings" button). Scale each trace separately or all traces with the same units together (in overlay mode: each trace scaled to fill the axis, or true values on one axis). RMS of what is on screen is shown under each name. Each trace has a "flip polarity" box that is remembered per channel and works together with whatever filter is selected (the global "Reverse polarity" flips everything). "Pick mode" turns clicks into event picks: click to add a pick, click on a pick line to remove it. Picks are saved automatically as UTC date-times in `miniseed_viewer_picks.csv` in the folder you started the program from (your home folder if that folder can't be written to; change the file in "Picks list..."); read them in Python with `pandas.read_csv(path, comment="#", parse_dates=["time_utc"])`. "Stats" adds a small plot under each trace of moving RMS, skewness and/or kurtosis (window length in seconds, default 5; not available in overlay mode). Zoom with the wheel or click twice on a plot; zooming out always shows all the data in the new interval. |
| Spectra | Power spectral density of the traces in the Magnify tab's current view. One spectrum for the whole view, or Welch average. |
| Spectrogram | Power (dB) versus time and frequency for each selected trace over the Magnify view (only the time axis zooms; frequency stays fixed). Window length (default 1024 samples), overlap, log/linear frequency, color map and range are adjustable; "Normalize each trace to its own peak" makes every trace's peak the top of the color scale. |
| Split | Spectra above, magnify below; spectra follow the magnify zoom. |

Mouse: left click = pick time, right drag = pan, wheel = zoom time.
"Settings" (top right) holds traces per screen, trace colors (color-blind-friendly by default, editable), channel sort order, time axis
(clock time or seconds since start of trace), default detrend, spectra defaults.
Settings, saved filters and calibrations are kept in `~/.config/miniseed_viewer/`
(`%APPDATA%\miniseed_viewer` on Windows, `~/Library/Application Support/miniseed_viewer` on Mac).

Notes on what the numbers mean:
* Everything is in raw counts (no instrument response).
* Spectra are one-sided power spectral density, counts^2/Hz, Hann taper, mean removed.
* A gap in the data breaks the line; spectra use the longest gap-free piece in view.
* Samples are held as 32-bit integers exactly as recorded. Filtered/detrended copies are
  made only for the traces on the current page and only for the chosen window.

## Installing

**Ubuntu** (22.04 or newer, 64-bit Intel/AMD): unpack `miniseed_viewer-VERSION-linux-x86_64.tar.gz` anywhere and run
`./miniseed_viewer file1.mseed ...` (or put the file somewhere on your PATH). It is one self-contained file;
it only needs the graphics library that every desktop Ubuntu already has (`libopengl0`).

**macOS** (11 or newer, Intel and Apple Silicon): unpack `miniseed_viewer-VERSION-macos.zip` and read
"READ ME FIRST.txt" inside (macOS asks for permission once because the app has no paid Apple signature).

Both files are produced by the automatic build in `.github/workflows/build.yml` (see "Making a release").

## Making a release

The files are built on GitHub's computers, so no Mac is needed:
1. Put the repository on GitHub (once): create an empty repository there, then
   `git remote add origin git@github.com:YOU/REPO.git` and `git push -u origin main`.
2. Every push builds Ubuntu and macOS versions (Actions tab > latest run > Artifacts).
3. To publish: change `VERSION`, commit, then `git tag v0.1.0 && git push origin v0.1.0`.
   The files appear under the repository's Releases page.

To build by hand: `packaging/build_linux.sh` (on the oldest Ubuntu you want to support) or
`packaging/build_macos.sh` (on a Mac). Results go to `dist/`.

## Building (Ubuntu)

    sudo apt install cmake g++ libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libgl-dev
    cmake -S . -B build && cmake --build build -j
    build/miniseed_viewer tests/data/*.mseed

All libraries are bundled as source in `third_party/` (see `third_party/README.md`),
so building needs no network access.

## Tests

    python tests/check_core.py build/test_core tests/data   # reading, filters, spectra, moving stats vs ObsPy/SciPy
    tests/run_ui.sh SCRIPT OUTDIR files...                          # drives the real window, saves screenshots

`run_ui.sh` plays back a text script (`frames N`, `click X Y`, `move`, `down`, `up`,
`rdown`, `rup`, `wheel DY`, `shot FILE.ppm`, `quit`) so the interface can be exercised
without a person. It uses a private settings folder, so it does not touch your own.

## Layout

    src/data.*      reading files, storage, fast min/max drawing summaries
    src/dsp.*       Butterworth filters, detrend, power spectra
    src/stats.*     moving RMS / skewness / kurtosis
    src/settings.*  saved options, saved filters, channel sort order
    src/app*.cpp    the interface (app_views: magnify/spectra/split; app_specgram: spectrogram; app_analyses: analyses tab)
    src/analysis.*  difference, cross-correlation and correlogram (no graphics; checked by `test_core analysis`)
    src/main.cpp    window setup (and the test-script player)
    packaging/      build scripts for the downloadable Ubuntu and macOS files
