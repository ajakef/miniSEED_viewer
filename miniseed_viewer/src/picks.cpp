#include "picks.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

#include "data.h"

namespace msv {

std::string format_iso(double t) {
    std::string s = format_time(t, 6);  // "YYYY-MM-DD HH:MM:SS.ffffff"
    s[10] = 'T';
    return s + "Z";
}

bool parse_iso(const std::string& in, double& t) {
    int y, mo, d, h, mi;
    double sec;
    std::string s = in;
    while (!s.empty() && (s.back() == 'Z' || s.back() == 'z' || s.back() == ' ' || s.back() == '\r')) s.pop_back();
    if (s.size() < 19) return false;
    if (s[10] == ' ') s[10] = 'T';
    if (std::sscanf(s.c_str(), "%d-%d-%dT%d:%d:%lf", &y, &mo, &d, &h, &mi, &sec) != 6) return false;
    if (mo < 1 || mo > 12 || d < 1 || d > 31 || h < 0 || h > 23 || mi < 0 || mi > 59 || sec < 0 || sec >= 61) return false;
    // Days since 1970-01-01 (Howard Hinnant's algorithm).
    y -= mo <= 2;
    const long long era = (y >= 0 ? y : y - 399) / 400;
    const long long yoe = y - era * 400;
    const long long doy = (153 * (mo + (mo > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    const long long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    const long long days = era * 146097 + doe - 719468;
    t = (double)days * 86400.0 + h * 3600.0 + mi * 60.0 + sec;
    return true;
}

bool save_picks(const std::string& path, std::vector<Pick> picks, std::string& error) {
    std::sort(picks.begin(), picks.end(), [](const Pick& a, const Pick& b) { return a.t < b.t; });
    std::ofstream f(path);
    if (!f) {
        error = "cannot write " + path;
        return false;
    }
    f << "# miniseed_viewer event picks. Times are UTC. In Python: "
         "pandas.read_csv(path, comment=\"#\", parse_dates=[\"time_utc\"])\n";
    f << "time_utc,channel,label\n";
    for (const Pick& p : picks) {
        std::string label = p.label;
        std::replace(label.begin(), label.end(), ',', ';');
        std::replace(label.begin(), label.end(), '\n', ' ');
        f << format_iso(p.t) << "," << p.id << "," << label << "\n";
    }
    f.flush();
    if (!f) {
        error = "error while writing " + path;
        return false;
    }
    return true;
}

bool load_picks(const std::string& path, std::vector<Pick>& picks, std::string& error) {
    std::ifstream f(path);
    if (!f) {
        error = "cannot read " + path;
        return false;
    }
    std::vector<Pick> out;
    std::string line;
    int bad = 0;
    while (std::getline(f, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
        if (line.empty() || line[0] == '#' || line.rfind("time_utc", 0) == 0) continue;
        size_t a = line.find(','), b = a == std::string::npos ? a : line.find(',', a + 1);
        Pick p;
        if (a == std::string::npos || b == std::string::npos || !parse_iso(line.substr(0, a), p.t)) {
            bad++;
            continue;
        }
        p.id = line.substr(a + 1, b - a - 1);
        p.label = line.substr(b + 1);
        out.push_back(p);
    }
    picks = out;
    error = bad ? std::to_string(bad) + " unreadable line(s) skipped" : "";
    return true;
}

}  // namespace msv
