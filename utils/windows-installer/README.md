<!--
SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
SPDX-License-Identifier: CC-BY-4.0
-->

# Windows Installer for Open-TYNDP

Guide for building and using the Windows installer based on pixi and NSIS.

**Table of Contents:**
- [Quick Start](#quick-start) - Get building quickly
- [Overview](#overview) - Understanding the approach
- [Building](#building-the-installer) - Detailed build instructions
- [Using](#using-the-installed-environment) - For end users
- [Customization](#customization) - Modify the installer
- [Troubleshooting](#troubleshooting) - Common issues and solutions
- [Advanced](#advanced-topics) - Cross-platform builds, automation, CI/CD

---

## Quick Start

**Automated Build (Recommended):**

The Windows installer is built automatically on every GitHub release. Simply create a new release tag and the installer will be attached as a release asset.

**Manual Build:**

```bash
# Navigate to installer directory
cd utils/windows-installer

# Run the build script
./build_pixi_installer.sh  # Linux/macOS/Git Bash

# Output: open-tyndp-0.4.0-pixi-Windows-x86_64.exe (~20 MB)
```

**File structure:**
```
utils/windows-installer/
├── build_pixi_installer.sh
├── pixi_installer.nsi
├── oet_logo.bmp
└── oet_logo.ico
```

---

## Overview

### What This Approach Does

Combines three technologies:
1. **pixi executable** - Bundled package manager (~20 MB)
2. **NSIS** - Creates Windows installer executable
3. **Git integration** - Automatically clones repository during installation

### Key Features

✅ **Small installer** - Only ~20 MB download (pixi executable only)
✅ **On-demand environment** - Downloads packages from conda-forge during installation
✅ **Integrated Git cloning** - One-click setup with custom directory selection
✅ **Cross-platform builds** - Build Windows installers from Linux/macOS
✅ **Direct shortcuts** - PowerShell and CMD shortcuts with embedded commands
✅ **No admin required** - User-level installation
✅ **Official branding** - Open Energy Transition logos and colors

### How It Works

1. **Installation**: Copies pixi.exe to `%LOCALAPPDATA%\open-tyndp`
2. **Repository Setup**: Clones repository using bundled pixi (via `pixi exec git`)
3. **Environment Installation**: Runs `pixi install` to download and set up conda environment
4. **Shortcuts**: Creates Start Menu shortcuts that activate environment via `pixi shell`

### Environment Activation

The shortcuts use `pixi shell` to activate the environment:
- **PowerShell**: `pixi shell -e open-tyndp` launches PowerShell in activated environment
- **CMD**: `pixi shell -e open-tyndp` launches CMD in activated environment
- **No PATH modification** - Environment only active in launched shell
- **No base environment** - Only the project environment is installed

---

## Building the Installer

### Automated Build

The installer is built automatically by GitHub Actions:
- **On every push/PR**: Workflow builds installer and uploads as artifact
- **On release tags**: Workflow builds installer and attaches to GitHub Release

### Using the Build Script

```bash
cd utils/windows-installer
./build_pixi_installer.sh  # Linux/macOS/Git Bash
```

The script automatically downloads pixi.exe if needed and builds the installer.

---

## How the Installer Works

### Installation Wizard Flow

1. **Welcome Page** - Introduction and license acceptance
2. **Repository Setup Page**
   - Checkbox to enable/disable auto-clone
   - Directory picker for repository location (default: `%USERPROFILE%\open-tyndp`)
   - Browse button for custom location
   - Validation of existing directories
   - Shows where pixi.exe will be installed
3. **Installation** - Progress bar showing:
   - Installing pixi executable
   - Cloning repository using `pixi exec git` (~1-2 minutes if enabled)
   - Installing environment using `pixi install` (~5-15 minutes, ~500-800 MB download)
4. **Finish** - Summary with next steps

### What Gets Installed

**Pixi executable** (`%LOCALAPPDATA%\open-tyndp`):
```
open-tyndp/
├── pixi.exe                  # Pixi package manager (~20 MB)
└── Uninstall.exe            # Uninstaller
```

**Repository files** (`%USERPROFILE%\open-tyndp` - if cloned):
```
open-tyndp/
├── .git/
├── .pixi/
│   └── envs/
│       └── open-tyndp/      # Conda environment installed by pixi
├── workflow/
├── config/
├── pixi.toml
├── pixi.lock
└── ... (project files)
```

**Start Menu shortcuts** (`%APPDATA%\Microsoft\Windows\Start Menu\Programs\Open-TYNDP`):
- Open-TYNDP PowerShell.lnk
- Open-TYNDP Command Prompt.lnk
- Uninstall Open-TYNDP.lnk

**Registry entries**:
- `HKCU\Software\Microsoft\Windows\CurrentVersion\Uninstall\Open-TYNDP` - For Add/Remove Programs
- `HKCU\Software\Open-TYNDP\RepositoryPath` - Stores repository location


---

## Using the Installed Environment

### For End Users

#### Installation

1. **Run the installer**
   - Double-click `open-tyndp-0.4.0-Windows-x86_64.exe`
   - Accept UAC prompt if needed

2. **Follow the wizard**
   - Click "Next" through Welcome and License
   - Choose environment location (or keep default)
   - **On Repository Setup page:**
     - ✅ Keep "Clone repository automatically" checked
     - Choose location or browse (default: `C:\Users\YourName\open-tyndp`)
     - Click "Next"
   - Wait for installation (3-5 minutes total)

3. **Start working**
   - Open Start Menu → "Open-TYNDP" → "Open-TYNDP PowerShell"
   - PowerShell opens in repository directory with environment activated
   - Run: `snakemake --help` to verify
   - Everything is ready!


### Uninstalling

**Via Add/Remove Programs:**
1. Open Settings → Apps → Installed apps
2. Find "Open-TYNDP"
3. Click "Uninstall"
4. **Important prompt:** "Remove repository directory?"
   - **Yes** - Deletes everything (environment + repository + your work!)
   - **No** - Keeps repository and your files, removes only environment

**Via Start Menu:**
- Start Menu → Open-TYNDP → Uninstall

---

## Troubleshooting


### Debugging Installation Issues

Enable NSIS logging (edit `.nsi`):

```nsis
# At top of file
!define ENABLE_LOGGING

# Installer will create install.log in installation directory
# View with: notepad %LOCALAPPDATA%\open-tyndp\install.log
```

---

## References and Resources

### Documentation
- [NSIS Documentation](https://nsis.sourceforge.io/Docs/) - Complete NSIS reference
- [NSIS Modern UI](https://nsis.sourceforge.io/Docs/Modern%20UI%202/Readme.html) - UI customization

### Examples and Templates
- [Conda Constructor NSIS Template](https://github.com/conda/constructor/blob/main/constructor/nsis/main.nsi.tmpl) - Reference implementation

### Tools
- [NSIS Download](https://nsis.sourceforge.io/Download) - Windows installer creator

### Support
- Project issues: https://github.com/open-energy-transition/open-tyndp/issues
- NSIS forums: https://forums.winamp.com/forum/151-nsis-discussion/
