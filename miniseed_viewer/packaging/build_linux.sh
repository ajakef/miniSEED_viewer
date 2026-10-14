#!/bin/bash
# Builds miniseed_viewer for Linux and packs it into dist/miniseed_viewer-VERSION-linux-x86_64.tar.gz
# Build this on the OLDEST Ubuntu you want to support (the program then also runs on newer ones).
set -e
cd "$(dirname "$0")/.."
version=$(cat VERSION)
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release -j"$(nproc)"
name="miniseed_viewer-$version-linux-$(uname -m)"
rm -rf "dist/$name"
mkdir -p "dist/$name"
cp build-release/miniseed_viewer README.md LICENSES.txt VERSION "dist/$name/"
strip "dist/$name/miniseed_viewer"
tar czf "dist/$name.tar.gz" -C dist "$name"
echo "built dist/$name.tar.gz (needs glibc >= $(objdump -T "dist/$name/miniseed_viewer" | grep -o 'GLIBC_[0-9.]*' | sort -uV | tail -1 | cut -d_ -f2))"
