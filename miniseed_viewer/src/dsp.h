// Filtering, detrending and power spectra.
#pragma once
#include <string>
#include <vector>

namespace msv {

enum class FilterType { LowPass = 0, HighPass = 1, BandPass = 2, BandStop = 3 };

// A Butterworth filter. f1 is the corner frequency for low/high pass; for band
// filters f1 and f2 are the lower and upper corners.
struct FilterSpec {
    std::string name = "new filter";
    FilterType type = FilterType::BandPass;
    int order = 4;
    double f1 = 1.0;
    double f2 = 10.0;
    bool zero_phase = true;  // run forward and backward (no time shift, double the order)
};

enum class Detrend { None = 0, Mean = 1, Linear = 2 };

// Empty string if the filter can be used at this sample rate, otherwise why not.
std::string check_filter(const FilterSpec& f, double rate);

// Filters x in place. Does nothing if check_filter() reports a problem.
void apply_filter(const FilterSpec& f, double rate, std::vector<double>& x);

void detrend(std::vector<double>& x, Detrend d);

struct Spectrum {
    std::vector<double> f;  // Hz (the zero-frequency point is left out)
    std::vector<double> p;  // power spectral density, (counts^2)/Hz
};

// Power spectral density of x. With welch = false the whole record gets one
// Hann-tapered FFT. With welch = true, windows of nwin samples (50% overlap)
// are averaged. Either way the mean is removed first.
Spectrum power_spectrum(const double* x, size_t n, double rate, bool welch, int nwin);

}  // namespace msv
