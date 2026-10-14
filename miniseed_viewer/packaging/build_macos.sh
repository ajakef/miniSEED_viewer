#!/bin/bash
# Builds miniseed_viewer for macOS (Intel + Apple Silicon in one file) and packs it into
# dist/miniseed_viewer-VERSION-macos.zip. No Apple developer account is needed: the app is
# signed "ad hoc" (free), which is enough for it to run once the user has allowed it.
set -e
cd "$(dirname "$0")/.."
version=$(cat VERSION)
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_ARCHITECTURES="arm64;x86_64"
cmake --build build-release -j"$(sysctl -n hw.ncpu)"

name="miniseed_viewer-$version-macos"
rm -rf "dist/$name"
app="dist/$name/miniseed_viewer.app"
mkdir -p "$app/Contents/MacOS"
cp build-release/miniseed_viewer "$app/Contents/MacOS/"
cat > "$app/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
  <key>CFBundleName</key><string>miniseed_viewer</string>
  <key>CFBundleDisplayName</key><string>miniseed_viewer</string>
  <key>CFBundleIdentifier</key><string>org.miniseed-viewer.app</string>
  <key>CFBundleExecutable</key><string>miniseed_viewer</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleVersion</key><string>$version</string>
  <key>CFBundleShortVersionString</key><string>$version</string>
  <key>LSMinimumSystemVersion</key><string>11.0</string>
  <key>NSHighResolutionCapable</key><true/>
</dict></plist>
PLIST
codesign --force --deep --sign - "$app"      # free "ad hoc" signature
cp README.md LICENSES.txt VERSION "dist/$name/"
cat > "dist/$name/READ ME FIRST.txt" <<TXT
miniseed_viewer $version for macOS (Intel and Apple Silicon)

This app is not signed with a paid Apple developer certificate, so macOS asks for permission once:
  1. Move miniseed_viewer.app to your Applications folder (or anywhere you like).
  2. Right-click (or Control-click) the app and choose Open, then press Open in the dialog.
     If macOS only offers "Move to Trash", open System Settings > Privacy & Security, scroll down
     and press "Open Anyway" next to the message about miniseed_viewer.
  Or, in Terminal, once:   xattr -dr com.apple.quarantine /path/to/miniseed_viewer.app

Use: drag miniSEED files onto the window, or from Terminal:
  /path/to/miniseed_viewer.app/Contents/MacOS/miniseed_viewer file1.mseed file2.mseed
(See README.md for how to use the program.)
TXT
rm -f "dist/$name.zip"
(cd dist && ditto -c -k --keepParent "$name" "$name.zip")
echo "built dist/$name.zip"
