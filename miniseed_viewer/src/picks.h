// Event picks: times marked by the user, saved in a plain text file.
#pragma once
#include <string>
#include <vector>

namespace msv {

struct Pick {
    double t = 0;       // seconds since 1970 (UTC)
    std::string id;     // channel the pick was made on, NET.STA.LOC.CHA
    std::string label;  // short note, e.g. "P"
};

// "2026-10-10T01:25:33.250000Z" (ISO 8601, UTC, microseconds).
std::string format_iso(double t);
// Reads that format (also without the fraction or the Z, or with a space instead of the T).
bool parse_iso(const std::string& s, double& t);

// File format: comment lines start with #, then a header line and one pick per line:
//   time_utc,channel,label
//   2026-10-10T01:25:33.250000Z,.10024..HDF,P
// In Python: pandas.read_csv(path, comment="#", parse_dates=["time_utc"])
bool save_picks(const std::string& path, std::vector<Pick> picks, std::string& error);
bool load_picks(const std::string& path, std::vector<Pick>& picks, std::string& error);

}  // namespace msv
