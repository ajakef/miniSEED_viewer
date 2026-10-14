#include "settings.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

namespace msv {

const char* sort_key_name(SortKey k) {
    switch (k) {
        case SortKey::CommandLine: return "Order files were listed";
        case SortKey::FileName: return "File name (alphabetical)";
        case SortKey::StartTime: return "Start time";
        case SortKey::Network: return "Network";
        case SortKey::Station: return "Station";
        case SortKey::Location: return "Location";
        case SortKey::Channel: return "Channel";
    }
    return "";
}

Config::Config() {
    // Default: files in the order given, then channel names as they were found.
    sort = {{SortKey::CommandLine, true, false}, {SortKey::FileName, false, false},
            {SortKey::StartTime, false, false},  {SortKey::Network, false, false},
            {SortKey::Station, false, false},    {SortKey::Location, false, false},
            {SortKey::Channel, false, false}};
}

std::string config_dir() {
    fs::path p;
#if defined(_WIN32)
    const char* base = std::getenv("APPDATA");
    p = base ? fs::path(base) : fs::path(".");
#elif defined(__APPLE__)
    const char* home = std::getenv("HOME");
    p = fs::path(home ? home : ".") / "Library" / "Application Support";
#else
    const char* xdg = std::getenv("XDG_CONFIG_HOME");
    const char* home = std::getenv("HOME");
    if (xdg && *xdg)
        p = xdg;
    else
        p = fs::path(home ? home : ".") / ".config";
#endif
    p /= "miniseed_viewer";
    std::error_code ec;
    fs::create_directories(p, ec);
    return p.string();
}

// ---------------------------------------------------------------------------
// Config file: one "name = value" per line.

void save_config(const Config& c) {
    std::ofstream f(config_dir() + "/config.txt");
    if (!f) return;
    f << "n_per_screen = " << c.n_per_screen << "\n"
      << "main_global_norm = " << c.main_global_norm << "\n"
      << "main_norm_visible = " << c.main_norm_visible << "\n"
      << "relative_time = " << c.relative_time << "\n"
      << "detrend = " << c.detrend << "\n"
      << "reverse_polarity = " << c.reverse_polarity << "\n"
      << "welch = " << c.welch << "\n"
      << "welch_window = " << c.welch_window << "\n"
      << "spec_log_x = " << c.spec_log_x << "\n"
      << "spec_log_y = " << c.spec_log_y << "\n"
      << "spec_normalize = " << c.spec_normalize << "\n"
      << "overlay_magnify = " << c.overlay_magnify << "\n"
      << "overlay_spectra = " << c.overlay_spectra << "\n"
      << "split_overlay_magnify = " << c.split_overlay_magnify << "\n"
      << "split_overlay_spectra = " << c.split_overlay_spectra << "\n"
      << "ui_scale = " << c.ui_scale << "\n";
    // key number, used (0/1), descending (0/1), most important first
    f << "sort =";
    for (const SortItem& s : c.sort) f << " " << (int)s.key << "," << s.use << "," << s.descending;
    f << "\n";
}

void load_config(Config& c) {
    std::ifstream f(config_dir() + "/config.txt");
    std::string line;
    while (std::getline(f, line)) {
        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        auto trim = [](std::string s) {
            size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
            return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
        };
        std::string k = trim(line.substr(0, eq)), v = trim(line.substr(eq + 1));
        auto b = [&] { return v == "1"; };
        try {
            if (k == "n_per_screen") c.n_per_screen = std::max(1, std::min(64, std::stoi(v)));
            else if (k == "main_global_norm") c.main_global_norm = b();
            else if (k == "main_norm_visible") c.main_norm_visible = b();
            else if (k == "relative_time") c.relative_time = b();
            else if (k == "detrend") c.detrend = std::max(0, std::min(2, std::stoi(v)));
            else if (k == "reverse_polarity") c.reverse_polarity = b();
            else if (k == "welch") c.welch = b();
            else if (k == "welch_window") c.welch_window = std::max(16, std::stoi(v));
            else if (k == "spec_log_x") c.spec_log_x = b();
            else if (k == "spec_log_y") c.spec_log_y = b();
            else if (k == "spec_normalize") c.spec_normalize = b();
            else if (k == "overlay_magnify") c.overlay_magnify = b();
            else if (k == "overlay_spectra") c.overlay_spectra = b();
            else if (k == "split_overlay_magnify") c.split_overlay_magnify = b();
            else if (k == "split_overlay_spectra") c.split_overlay_spectra = b();
            else if (k == "ui_scale") c.ui_scale = std::max(0.5f, std::min(4.0f, std::stof(v)));
            else if (k == "sort") {
                std::vector<SortItem> items;
                std::istringstream ss(v);
                std::string tok;
                while (ss >> tok) {
                    int key, use, desc;
                    if (std::sscanf(tok.c_str(), "%d,%d,%d", &key, &use, &desc) == 3 && key >= 0 && key < kNumSortKeys)
                        items.push_back({(SortKey)key, use != 0, desc != 0});
                }
                // Accept only a complete list (every key exactly once).
                std::vector<bool> seen(kNumSortKeys, false);
                bool ok = (int)items.size() == kNumSortKeys;
                for (auto& it : items) {
                    if (seen[(int)it.key]) ok = false;
                    seen[(int)it.key] = true;
                }
                if (ok) c.sort = items;
            }
        } catch (...) {
        }
    }
}

// ---------------------------------------------------------------------------
// Filter file: one filter per line, "name|type|order|f1|f2|zero_phase".

static std::vector<FilterSpec> default_filters() {
    auto mk = [](const char* name, FilterType t, int order, double f1, double f2) {
        FilterSpec f;
        f.name = name;
        f.type = t;
        f.order = order;
        f.f1 = f1;
        f.f2 = f2;
        return f;
    };
    return {mk("High-pass 1 Hz", FilterType::HighPass, 4, 1, 0), mk("Low-pass 10 Hz", FilterType::LowPass, 4, 10, 0),
            mk("Band-pass 1-10 Hz", FilterType::BandPass, 4, 1, 10)};
}

void save_filters(const std::vector<FilterSpec>& filters) {
    std::ofstream f(config_dir() + "/filters.txt");
    if (!f) return;
    for (const FilterSpec& s : filters) {
        std::string name = s.name;
        std::replace(name.begin(), name.end(), '|', '/');
        std::replace(name.begin(), name.end(), '\n', ' ');
        f << name << "|" << (int)s.type << "|" << s.order << "|" << s.f1 << "|" << s.f2 << "|" << s.zero_phase << "\n";
    }
}

void load_filters(std::vector<FilterSpec>& filters) {
    filters.clear();
    std::ifstream f(config_dir() + "/filters.txt");
    if (!f) {
        filters = default_filters();
        save_filters(filters);
        return;
    }
    std::string line;
    while (std::getline(f, line)) {
        std::vector<std::string> parts;
        std::stringstream ss(line);
        std::string p;
        while (std::getline(ss, p, '|')) parts.push_back(p);
        if (parts.size() != 6) continue;
        try {
            FilterSpec s;
            s.name = parts[0];
            int t = std::stoi(parts[1]);
            if (t < 0 || t > 3) continue;
            s.type = (FilterType)t;
            s.order = std::max(1, std::min(10, std::stoi(parts[2])));
            s.f1 = std::stod(parts[3]);
            s.f2 = std::stod(parts[4]);
            s.zero_phase = parts[5] == "1";
            filters.push_back(s);
        } catch (...) {
        }
    }
}

// ---------------------------------------------------------------------------
// Calibration file: one channel per line, "id|scale|units".

void save_calib(const std::map<std::string, Calib>& calib) {
    std::ofstream f(config_dir() + "/calibration.txt");
    if (!f) return;
    f.precision(12);
    for (const auto& kv : calib) {
        if (kv.second.scale == 1.0 && kv.second.units == "counts") continue;  // nothing to remember
        std::string units = kv.second.units;
        std::replace(units.begin(), units.end(), '|', '/');
        f << kv.first << "|" << kv.second.scale << "|" << units << "\n";
    }
}

void load_calib(std::map<std::string, Calib>& calib) {
    calib.clear();
    std::ifstream f(config_dir() + "/calibration.txt");
    std::string line;
    while (std::getline(f, line)) {
        std::vector<std::string> parts;
        std::stringstream ss(line);
        std::string p;
        while (std::getline(ss, p, '|')) parts.push_back(p);
        if (parts.size() < 2) continue;
        try {
            Calib c;
            c.scale = std::stod(parts[1]);
            c.units = parts.size() > 2 && !parts[2].empty() ? parts[2] : "counts";
            if (c.scale != 0) calib[parts[0]] = c;
        } catch (...) {
        }
    }
}

// ---------------------------------------------------------------------------

std::vector<int> sort_order(const std::vector<Trace>& traces, const std::vector<SortItem>& order) {
    auto base = [](const std::string& p) {
        size_t i = p.find_last_of("/\\");
        return i == std::string::npos ? p : p.substr(i + 1);
    };
    std::vector<int> idx(traces.size());
    for (size_t i = 0; i < idx.size(); i++) idx[i] = (int)i;
    std::stable_sort(idx.begin(), idx.end(), [&](int ia, int ib) {
        const Trace& a = traces[ia];
        const Trace& b = traces[ib];
        for (const SortItem& s : order) {
            if (!s.use) continue;
            int c = 0;  // <0: a first
            switch (s.key) {
                case SortKey::CommandLine: c = (a.file_order > b.file_order) - (a.file_order < b.file_order); break;
                case SortKey::FileName: c = base(a.file).compare(base(b.file)); break;
                case SortKey::StartTime: c = (a.t_min > b.t_min) - (a.t_min < b.t_min); break;
                case SortKey::Network: c = a.net.compare(b.net); break;
                case SortKey::Station: c = a.sta.compare(b.sta); break;
                case SortKey::Location: c = a.loc.compare(b.loc); break;
                case SortKey::Channel: c = a.cha.compare(b.cha); break;
            }
            if (c != 0) return s.descending ? c > 0 : c < 0;
        }
        return false;
    });
    return idx;
}

}  // namespace msv
