#!/bin/bash

# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

set -e

# Get script directory and change to it
cd $(dirname "$0")

echo "Building Open-TYNDP Windows Installer ..."
echo ""

# Get version from git tag if not provided
if [ -z "$VERSION" ]; then
    VERSION=$(git describe --tags --always --dirty 2>/dev/null || echo "0.0.0-dev")
    VERSION=${VERSION#v}  # Remove 'v' prefix if present
fi
echo "Building installer version: $VERSION"
echo ""

# Create shallow clone for bundling
echo "Step 1/3: Creating shallow git clone for bundling..."
REPO_ROOT=$(git rev-parse --show-toplevel)

# Remove left-over repo-bundle
rm -fr repo-bundle

# Create shallow clone from current HEAD
git clone --depth 1 "file://$REPO_ROOT/.git" "repo-bundle"
(
    cd "repo-bundle"
    # Set remote to GitHub URL and expand the fetch refspec again
    # (so that after git fetch, git checkout/switch works normally)
    git remote set-url origin "${PRODUCT_REPO_URL:-https://github.com/open-energy-transition/open-tyndp.git}"
    git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"

    # Windows gets confused about file mode (executable bit or no), so let's disable that
    git config core.fileMode false
)

REPO_SIZE=$(du -sk repo-bundle | cut -f1)
echo "Created repo-bundle directory ($(du -sh repo-bundle | cut -f1))"
echo ""

# Download pixi executable for Windows if not present
if [ ! -f pixi.exe ]; then
    echo "Step 2/3: Downloading pixi executable for Windows..."
    curl "https://github.com/prefix-dev/pixi/releases/latest/download/pixi-x86_64-pc-windows-msvc.exe" -Lo pixi.exe
else
    echo "Step 2/3: Using existing pixi.exe"
fi
PIXI_SIZE=$(du -sk pixi.exe | cut -f1)
echo ""

# Calculate total estimated size (pixi.exe + repo + uninstaller overhead ~100KB)
ESTIMATED_SIZE=$((PIXI_SIZE + REPO_SIZE + 100))
echo "Estimated installation size: $((ESTIMATED_SIZE / 1024)) MB"
echo ""

# Build installer
echo "Step 3/3: Building installer..."
pixi exec --spec nsis makensis -DPRODUCT_VERSION=$VERSION -DESTIMATED_SIZE=$ESTIMATED_SIZE installer.nsi

echo ""
echo "✅ Build complete!"
echo "   Installer: open-tyndp-${VERSION}-Windows-x86_64.exe"

# Clean up bundle directory
rm -rf repo-bundle
