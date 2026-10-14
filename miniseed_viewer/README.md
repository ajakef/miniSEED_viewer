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
| Magnify | Filtering, detrending and polarity reversal (for all traces, or per trace with each trace's "settings" button). Scale each trace separately or all traces with the same units together. RMS of what is on screen is shown under each name. "Stats" adds a small plot under each trace of moving RMS, skewness and/or kurtosis (window length in seconds, default 5; not available in overlay mode). Zoom with the wheel or click twice on a plot; zooming out always shows all the data in the new interval. |
| Spectra | Power spectral density of the traces in the Magnify tab's current view. One spectrum for the whole view, or Welch average. |
| Spectrogram | Power (dB) versus time and frequency for each selected trace over the Magnify view (only the time axis zooms; frequency stays fixed). Window length (default 1024 samples), overlap, log/linear frequency, color map and range are adjustable. |
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
    src/app*.cpp    the interface (app_views: magnify/spectra/split; app_specgram: spectrogram)
    src/main.cpp    window setup (and the test-script player)
