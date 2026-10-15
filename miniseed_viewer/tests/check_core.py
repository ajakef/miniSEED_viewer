"""Compare the C++ core against ObsPy / SciPy.

usage: python check_core.py PATH_TO_test_core SAMPLE_DIR
"""
import glob
import re
import subprocess
import sys
import tempfile
import os

import numpy as np
import obspy
from scipy import signal

exe, sample_dir = sys.argv[1], sys.argv[2]
ok = True


def report(name, good, detail=""):
    global ok
    ok &= bool(good)
    print(("PASS " if good else "FAIL ") + name + ("  " + detail if detail else ""))


# ---- reading ----
files = sorted(glob.glob(os.path.join(sample_dir, "*.mseed")))
out = subprocess.run([exe, "load"] + files, capture_output=True, text=True).stdout
traces = {}
cur = None
for line in out.splitlines():
    m = re.match(r"TRACE (\S+) nseg=(\d+) rate=([\d.]+) t_min=([\d.]+)", line)
    if m:
        cur = m.group(1)
        traces[cur] = {"nseg": int(m.group(2)), "rate": float(m.group(3)), "t0": float(m.group(4)), "segs": []}
    m = re.match(r"\s+SEG t0=([\d.]+) n=(\d+) first=(\S+),(\S+),(\S+) last=(\S+) rms=(\S+)", line)
    if m:
        traces[cur]["segs"].append(m.groups())
    if "MINMAX_BAD" in line:
        report("min/max summary matches plain scan", line.split()[-1] == "0")
report("number of traces", len(traces) == len(files), f"{len(traces)} vs {len(files)}")
for f in files:
    st = obspy.read(f)
    tr = st[0]
    tid = f"{tr.stats.network}.{tr.stats.station}.{tr.stats.location}.{tr.stats.channel}"
    t = traces.get(tid)
    if t is None:
        report("trace " + tid, False, "missing; have " + ",".join(traces))
        continue
    seg = t["segs"][0]
    same = (
        int(seg[1]) == tr.stats.npts
        and abs(float(seg[0]) - tr.stats.starttime.timestamp) < 1e-5
        and [float(v) for v in seg[2:5]] == list(tr.data[:3].astype(float))
        and float(seg[5]) == float(tr.data[-1])
        and abs(float(seg[6]) - np.sqrt(np.mean(tr.data.astype(float) ** 2))) < 1e-3 * abs(float(seg[6]))
    )
    report("read " + tid, same and t["nseg"] == 1 and abs(t["rate"] - tr.stats.sampling_rate) < 1e-6)

# ---- filters and spectra ----
d = tempfile.mkdtemp()
txt = subprocess.run([exe, "dump", d], capture_output=True, text=True).stdout
print(txt.strip())
rd = lambda n: np.fromfile(os.path.join(d, n))
x = rd("x.bin")
fs = 200
mid = slice(len(x) // 10, -len(x) // 10)
cases = [("lowpass", 4, 10, None, True), ("highpass", 3, 1, None, True), ("bandpass", 4, 1, 10, True),
         ("bandstop", 2, 4, 6, True), ("bandpass", 3, 2, 20, False), ("lowpass", 1, 20, None, False),
         ("bandstop", 3, 4, 40, True), ("lowpass", 5, 3, None, True), ("bandpass", 5, 0.5, 30, True),
         ("highpass", 4, 0.1, None, False)]
for i, (btype, order, f1, f2, zp) in enumerate(cases):
    wn = f1 if f2 is None else [f1, f2]
    sos = signal.butter(order, wn, btype=btype, fs=fs, output="sos")
    padlen = min(max(int(np.ceil(3 * fs / f1)), 60), len(x) - 1)  # same padding as the C++ code
    ref = signal.sosfiltfilt(sos, x, padlen=padlen) if zp else signal.sosfilt(sos, x)
    mine = rd(f"y{i}.bin")
    sl = mid if zp else slice(len(x) // 2, None)  # causal filter start-up differs at the start
    err = np.max(np.abs(mine[sl] - ref[sl])) / np.max(np.abs(ref[sl]))
    report(f"filter {btype} order {order} zero_phase={zp}", err < 1e-6, f"relative error {err:.2e}")
yl = rd("detrend_linear.bin")
t = np.arange(len(x))
ref = signal.detrend(x, type="linear")
report("linear detrend", np.max(np.abs(yl - ref)) < 1e-6)
w = signal.windows.hann(len(x), sym=True)
f, p = signal.welch(x, fs=fs, window=w, nperseg=len(x), noverlap=0, detrend="constant", scaling="density")
fm, pm = rd("spec1_f.bin"), rd("spec1_p.bin")
report("single spectrum", np.allclose(fm, f[1:]) and np.allclose(pm, p[1:], rtol=1e-9), "")
w2 = signal.windows.hann(1024, sym=True)
f, p = signal.welch(x, fs=fs, window=w2, nperseg=1024, noverlap=512, detrend="constant", scaling="density")
fm, pm = rd("spec2_f.bin"), rd("spec2_p.bin")
report("Welch spectrum", np.allclose(fm, f[1:]) and np.allclose(pm, p[1:], rtol=1e-9), "")
# ---- moving statistics (integer data, 5 s window)
from scipy import stats as sps
xi = np.round(x)
t = rd("ms_t.bin")
centres = np.round((t - 100.0) * fs).astype(int)
W = 1000
rms_ref, sk_ref, ku_ref = [], [], []
for c in centres:
    w = xi[c - W // 2: c - W // 2 + W]
    rms_ref.append(np.sqrt(np.mean(w ** 2)))
    sk_ref.append(sps.skew(w))
    ku_ref.append(sps.kurtosis(w, fisher=False))
report("moving stats: some windows", len(t) > 100, f"{len(t)} windows")
report("moving RMS", np.allclose(rd("ms_rms.bin"), rms_ref, rtol=1e-9))
report("moving skewness", np.allclose(rd("ms_skew.bin"), sk_ref, rtol=1e-6, atol=1e-8))
report("moving kurtosis", np.allclose(rd("ms_kurt.bin"), ku_ref, rtol=1e-6, atol=1e-8))
# ---- picks file: readable in Python, ISO times
import csv
rows = [r for r in csv.reader(l for l in open(os.path.join(d, "picks.csv")) if not l.startswith("#"))]
report("picks file header", rows[0] == ["time_utc", "channel", "label"])
times = [obspy.UTCDateTime(r[0]).timestamp for r in rows[1:]]
report("picks sorted by time and in ISO format", times == sorted(times) and rows[1][0] == "2026-10-10T00:00:00.000000Z", rows[1][0])
report("pick times exact", abs(times[2] - 1791600000.123456) < 1e-5 and abs(times[1] - 1791597933.25) < 1e-6)
report("comma in label replaced", rows[1][2] == "a;b")
try:
    import pandas as pd
    df = pd.read_csv(os.path.join(d, "picks.csv"), comment="#", parse_dates=["time_utc"])
    report("pandas reads the picks file", len(df) == 3 and str(df.time_utc.dtype).startswith("datetime64"), str(df.time_utc.dtype))
except ImportError:
    print("(pandas not installed: skipped)")
print("ALL PASS" if ok else "SOME FAILED")
sys.exit(0 if ok else 1)
