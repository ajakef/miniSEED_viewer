#!/bin/bash
# Plays a script of mouse actions against the real window and saves screenshots (dev tool).
#   run_ui.sh SCRIPT.txt OUTDIR [files...]   (settings go to OUTDIR/home, not your real settings)
set -e
script=$(realpath "$1"); out=$(realpath -m "$2"); shift 2
mkdir -p "$out/home"
here=$(dirname "$(realpath "$0")")
HOME="$out/home" XDG_CONFIG_HOME="$out/home/.config" MSV_SCRIPT="$script" timeout 120 "$here/../build/miniseed_viewer" "$@"
for f in "$out"/*.ppm; do [ -e "$f" ] && convert "$f" "${f%.ppm}.png" && rm "$f"; done
