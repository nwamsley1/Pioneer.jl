#!/usr/bin/env bash
set -euo pipefail

# Download and extract the latest Poppler release for Windows
release_json=$(curl -s https://api.github.com/repos/oschwartz10612/poppler-windows/releases/latest)
asset_url=$(echo "$release_json" | grep browser_download_url | grep x86_64.zip | head -n 1 | cut -d '"' -f 4)

workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT

curl -Lf --retry 3 --retry-delay 2 "$asset_url" -o "$workdir/poppler.zip"

if command -v unzip >/dev/null 2>&1; then
    unzip -q "$workdir/poppler.zip" -d "$workdir"
elif command -v 7z >/dev/null 2>&1; then
    7z x "$workdir/poppler.zip" -o"$workdir" >/dev/null
else
    echo "No unzip utility found" >&2
    exit 1
fi

bindir=$(find "$workdir" -path '*pdfunite.exe' -exec dirname {} \; | head -n 1)

# Add Poppler bin directory to PATH for subsequent steps
if [ -n "$bindir" ]; then
    echo "$bindir" >> "$GITHUB_PATH"
fi