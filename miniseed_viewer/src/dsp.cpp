#include "dsp.h"

#include <algorithm>
#include <cmath>
#include <complex>

#include "pocketfft_hdronly.h"

namespace msv {

namespace {

using cd = std::complex<double>;
const double kPi = 3.14159265358979323846;

// One second-order (or first-order, with b2 = a2 = 0) filter section.
struct Sos {
    double b0, b1, b2, a1, a2;
};

double section_gain(const Sos& s, double theta) {
    cd z1 = std::polar(1.0, -theta), z2 = std::polar(1.0, -2 * theta);
    cd num = s.b0 + s.b1 * z1 + s.b2 * z2;
    cd den = 1.0 + s.a1 * z1 + s.a2 * z2;
    return std::abs(num / den);
}

// Butterworth design: analog prototype -> frequency transform -> bilinear
// transform -> pairs of conjugate poles as second-order sections.
std::vector<Sos> design(const FilterSpec& f, double fs) {
    const int N = f.order;
    const double fs2 = 2 * fs;
    auto warp = [&](double hz) { return fs2 * std::tan(kPi * hz / fs); };
    auto bilinear = [&](cd s) { return (fs2 + s) / (fs2 - s); };

    // Prototype poles, one from each conjugate pair (and the real pole if N is odd).
    std::vector<cd> proto;
    for (int k = 0; k < N; k++) {
        double th = kPi * (2.0 * k + N + 1) / (2.0 * N);
        cd p(std::cos(th), std::sin(th));
        if (p.imag() > 1e-9) proto.push_back(p);
    }
    const bool odd = (N % 2) == 1;

    std::vector<Sos> out;
    double ref_theta = 0;  // frequency at which every section gets gain 1
    double theta0 = 0;     // centre of a band filter (radians per sample)
    const bool band = f.type == FilterType::BandPass || f.type == FilterType::BandStop;
    double w0 = 0, bw = 0;
    if (band) {
        double W1 = warp(f.f1), W2 = warp(f.f2);
        w0 = std::sqrt(W1 * W2);
        bw = W2 - W1;
        theta0 = 2 * std::atan(w0 / fs2);
    }
    double Wc = band ? 0 : warp(f.f1);

    auto numerator = [&](bool second_order, double& b0, double& b1, double& b2) {
        switch (f.type) {
            case FilterType::LowPass:
                b0 = 1; b1 = second_order ? 2 : 1; b2 = second_order ? 1 : 0;
                break;
            case FilterType::HighPass:
                b0 = 1; b1 = second_order ? -2 : -1; b2 = second_order ? 1 : 0;
                break;
            case FilterType::BandPass:
                b0 = 1; b1 = 0; b2 = -1;
                break;
            case FilterType::BandStop:
                b0 = 1; b1 = -2 * std::cos(theta0); b2 = 1;
                break;
        }
    };

    auto add_pair = [&](cd zpole) {  // zpole and its conjugate
        Sos s;
        s.a1 = -2 * zpole.real();
        s.a2 = std::norm(zpole);
        numerator(true, s.b0, s.b1, s.b2);
        out.push_back(s);
    };
    auto add_real = [&](double r) {
        Sos s;
        s.a1 = -r;
        s.a2 = 0;
        numerator(false, s.b0, s.b1, s.b2);
        out.push_back(s);
    };

    if (!band) {
        auto map = [&](cd p) { return f.type == FilterType::LowPass ? Wc * p : Wc / p; };
        for (cd p : proto) add_pair(bilinear(map(p)));
        if (odd) add_real(bilinear(map(cd(-1, 0))).real());
        ref_theta = (f.type == FilterType::LowPass) ? 0 : kPi;
    } else {
        auto lp2b = [&](cd p, cd& q1, cd& q2) {
            cd a = (f.type == FilterType::BandPass) ? p * bw / 2.0 : (bw / 2.0) / p;
            cd d = std::sqrt(a * a - w0 * w0);
            q1 = a + d;
            q2 = a - d;
        };
        for (cd p : proto) {
            cd q1, q2;
            lp2b(p, q1, q2);
            add_pair(bilinear(q1));
            add_pair(bilinear(q2));
        }
        if (odd) {
            // The real prototype pole becomes two poles that are either a
            // conjugate pair or (for a wide band) two real poles.
            double a = (f.type == FilterType::BandPass) ? -bw / 2.0 : -(bw / 2.0);
            double disc = a * a - w0 * w0;
            if (disc >= 0) {
                double d = std::sqrt(disc);
                double r1 = bilinear(cd(a + d, 0)).real(), r2 = bilinear(cd(a - d, 0)).real();
                Sos s;
                s.a1 = -(r1 + r2);
                s.a2 = r1 * r2;
                numerator(true, s.b0, s.b1, s.b2);
                out.push_back(s);
            } else {
                add_pair(bilinear(cd(a, std::sqrt(-disc))));
            }
        }
        ref_theta = (f.type == FilterType::BandPass) ? theta0 : 0;
    }

    for (Sos& s : out) {
        double g = section_gain(s, ref_theta);
        if (g > 0) {
            s.b0 /= g;
            s.b1 /= g;
            s.b2 /= g;
        }
    }
    return out;
}

void run_sos(const std::vector<Sos>& sections, double* x, size_t n) {
    for (const Sos& s : sections) {
        double z1 = 0, z2 = 0;
        for (size_t i = 0; i < n; i++) {
            double xi = x[i];
            double y = s.b0 * xi + z1;
            z1 = s.b1 * xi - s.a1 * y + z2;
            z2 = s.b2 * xi - s.a2 * y;
            x[i] = y;
        }
    }
}

}  // namespace

std::string check_filter(const FilterSpec& f, double rate) {
    if (rate <= 0) return "unknown sample rate";
    if (f.order < 1 || f.order > 10) return "order must be 1 to 10";
    double nyq = rate / 2;
    if (!(f.f1 > 0) || f.f1 >= nyq) return "corner frequency must be between 0 and " + std::to_string(nyq) + " Hz";
    if (f.type == FilterType::BandPass || f.type == FilterType::BandStop) {
        if (!(f.f2 > f.f1)) return "upper corner must be above lower corner";
        if (f.f2 >= nyq) return "upper corner must be below " + std::to_string(nyq) + " Hz";
    }
    return "";
}

void apply_filter(const FilterSpec& f, double rate, std::vector<double>& x) {
    if (!check_filter(f, rate).empty() || x.size() < 4) return;
    std::vector<Sos> sos = design(f, rate);
    if (!f.zero_phase) {
        run_sos(sos, x.data(), x.size());
        return;
    }
    // Forward-backward filtering. Pad both ends by odd reflection so the
    // start-up wiggle falls in the padding, which is thrown away.
    const size_t n = x.size();
    size_t pad = (size_t)std::ceil(3.0 * rate / f.f1);
    pad = std::max<size_t>(pad, 60);
    pad = std::min(pad, n - 1);
    std::vector<double> w(n + 2 * pad);
    for (size_t i = 0; i < pad; i++) {
        w[pad - 1 - i] = 2 * x[0] - x[i + 1];
        w[pad + n + i] = 2 * x[n - 1] - x[n - 2 - i];
    }
    std::copy(x.begin(), x.end(), w.begin() + pad);
    run_sos(sos, w.data(), w.size());
    std::reverse(w.begin(), w.end());
    run_sos(sos, w.data(), w.size());
    std::reverse(w.begin(), w.end());
    std::copy(w.begin() + pad, w.begin() + pad + n, x.begin());
}

void detrend(std::vector<double>& x, Detrend d) {
    const size_t n = x.size();
    if (d == Detrend::None || n == 0) return;
    double mean = 0;
    for (double v : x) mean += v;
    mean /= (double)n;
    if (d == Detrend::Mean || n < 2) {
        for (double& v : x) v -= mean;
        // Second pass removes the rounding error of the first, so constant data become exactly zero.
        double resid = 0;
        for (double v : x) resid += v;
        resid /= (double)n;
        for (double& v : x) v -= resid;
        return;
    }
    // Least-squares line, with the index centred for accuracy.
    double c = 0.5 * (double)(n - 1), sxy = 0, sxx = 0;
    for (size_t i = 0; i < n; i++) {
        double dx = (double)i - c;
        sxy += dx * (x[i] - mean);
        sxx += dx * dx;
    }
    double slope = sxx > 0 ? sxy / sxx : 0;
    for (size_t i = 0; i < n; i++) x[i] -= mean + slope * ((double)i - c);
}

namespace {

// Adds the power of one demeaned, Hann-tapered block to acc (length m/2+1).
void add_block_power(const double* x, size_t m, std::vector<double>& acc, double& wsum2) {
    std::vector<double> buf(m);
    double mean = 0;
    for (size_t i = 0; i < m; i++) mean += x[i];
    mean /= (double)m;
    wsum2 = 0;
    for (size_t i = 0; i < m; i++) {
        double w = m > 1 ? 0.5 - 0.5 * std::cos(2 * kPi * (double)i / (double)(m - 1)) : 1.0;
        buf[i] = (x[i] - mean) * w;
        wsum2 += w * w;
    }
    std::vector<cd> X(m / 2 + 1);
    pocketfft::r2c(pocketfft::shape_t{m}, pocketfft::stride_t{(ptrdiff_t)sizeof(double)},
                   pocketfft::stride_t{(ptrdiff_t)sizeof(cd)}, 0, true, buf.data(), X.data(), 1.0);
    for (size_t k = 0; k < X.size(); k++) acc[k] += std::norm(X[k]);
}

}  // namespace

Spectrum power_spectrum(const double* x, size_t n, double rate, bool welch, int nwin) {
    Spectrum s;
    if (n < 8 || rate <= 0) return s;
    size_t m = n;
    size_t nseg = 1, hop = n;
    if (welch) {
        size_t w = (size_t)std::max(nwin, 16);
        if (w < n) {
            m = w;
            hop = std::max<size_t>(m / 2, 1);
            nseg = (n - m) / hop + 1;
        }
    }
    std::vector<double> acc(m / 2 + 1, 0.0);
    double wsum2 = 1;
    for (size_t k = 0; k < nseg; k++) add_block_power(x + k * hop, m, acc, wsum2);
    const double scale = 1.0 / ((double)nseg * rate * wsum2);
    const size_t nk = acc.size();
    s.f.reserve(nk);
    s.p.reserve(nk);
    for (size_t k = 1; k < nk; k++) {
        double p = acc[k] * scale;
        bool nyquist = (m % 2 == 0) && (k == nk - 1);
        if (!nyquist) p *= 2;  // one-sided
        s.f.push_back((double)k * rate / (double)m);
        s.p.push_back(p);
    }
    return s;
}

}  // namespace msv
