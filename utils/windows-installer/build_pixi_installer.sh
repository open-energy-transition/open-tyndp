#!/bin/bash

# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

set -e

echo "Building Open-TYNDP Windows Installer ..."
echo ""

# Get version from git tag if not provided
if [ -z "$VERSION" ]; then
    VERSION=$(git describe --tags --always --dirty 2>/dev/null || echo "0.0.0-dev")
    VERSION=${VERSION#v}  # Remove 'v' prefix if present
fi
echo "Building installer version: $VERSION"
echo ""

# Download pixi executable for Windows if not present
if [ ! -f pixi.exe ]; then
    echo "Step 1/2: Downloading pixi executable for Windows..."
    PIXI_VERSION=$(curl -s https://api.github.com/repos/prefix-dev/pixi/releases/latest | grep '"tag_name"' | sed -E 's/.*"v([^"]+)".*/\1/')
    echo "Latest pixi version: v$PIXI_VERSION"
    wget "https://github.com/prefix-dev/pixi/releases/latest/download/pixi-x86_64-pc-windows-msvc.exe" -O pixi.exe
    echo "Downloaded pixi.exe ($(du -h pixi.exe | cut -f1))"
else
    echo "Step 1/2: Using existing pixi.exe ($(du -h pixi.exe | cut -f1))"
fi

# Build installer
echo "Step 2/2: Building installer..."
pixi exec --spec nsis makensis -DPRODUCT_VERSION=$VERSION pixi_installer.nsi

echo ""
echo "✅ Build complete!"
echo "   Installer: open-tyndp-${VERSION}-pixi-Windows-x86_64.exe"
ls -lh open-tyndp-*-pixi-Windows-x86_64.exe
