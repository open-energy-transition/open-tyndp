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
./build_installer.sh  # Linux/macOS/Git Bash

# Output: open-tyndp-0.4.0-pixi-Windows-x86_64.exe (~20 MB)
```

**File structure:**
```
utils/windows-installer/
├── build_installer.sh
├── installer.nsi
├── oet_logo.bmp
└── oet_logo.ico
```

---

## Overview

### Key Features

- Bundles pixi executable and git repository in a single installer using NSIS
- Downloads conda packages during installation (~500-800 MB)
- Can be built on Linux/macOS for Windows targets
- Includes PowerShell and Command Prompt shortcuts
- Requires no administrator rights
- Uses Open Energy Transition branding

### How It Works

1. Copies pixi.exe to `%LOCALAPPDATA%\open-tyndp`
2. Extracts bundled repository archive to chosen location
3. Runs `pixi install` to download and set up conda environment
4. Creates Start Menu shortcuts that activate environment via `pixi shell`

### Environment Activation

The shortcuts use `pixi shell` to activate the environment:
- Launches PowerShell or Command Prompt with environment activated
- No system PATH modification, environment only active in launched shell

---

## Building the Installer

### Automated Build (Recommended)

The installer is built automatically by GitHub Actions on release tags and attached to GitHub Releases.

### Manual Build

```bash
cd utils/windows-installer
./build_installer.sh  # Linux/macOS/Git Bash
```

The script automatically:
1. Downloads pixi.exe if needed
2. Creates a git archive of the repository
3. Builds the NSIS installer (~50-60 MB)

---

## How the Installer Works

### Installation Wizard Flow

1. Welcome page with introduction and license acceptance
2. Repository setup page with directory picker and path validation (no spaces allowed)
3. Installation progress showing extraction and environment setup (~5-15 minutes)
4. Finish page with summary

### What Gets Installed

**Pixi executable** (`%LOCALAPPDATA%\open-tyndp`):
```
open-tyndp/
├── pixi.exe                  # Pixi package manager (~60 MB)
└── Uninstall.exe            # Uninstaller
```

**Repository files** (`%USERPROFILE%\open-tyndp` - if extracted):
```
open-tyndp/
├── .git/                    # Full git repository
├── .pixi/
│   └── envs/
│       └── open-tyndp/      # Conda environment (~3 GB)
├── rules/
├── scripts/
├── config/
├── pixi.toml
└── ... (all project files)
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
   - Double-click `open-tyndp-0.4.0-pixi-Windows-x86_64.exe`
   - No admin rights required

2. **Follow the wizard**
   - Click "Next" through Welcome and License
   - **On Repository Setup page:**
     - ✅ Keep "Extract repository and install environment automatically" checked
     - Choose location without spaces (default: `C:\Users\YourName\open-tyndp`)
     - Click "Next"
   - Wait for installation (~5-15 minutes depending on internet speed)

3. **Start working**
   - Open Start Menu → "Open-TYNDP" → "Open-TYNDP PowerShell"
   - PowerShell opens in repository directory with environment activated
   - Run: `snakemake --help` to verify


### Uninstalling

**Via Add/Remove Programs or Start Menu:**
1. Settings → Apps → Installed apps → "Open-TYNDP" → Uninstall
   - OR: Start Menu → Open-TYNDP → Uninstall Open-TYNDP
2. **Choose what to remove:**
   - ☐ Remove repository directory (⚠️ deletes all your work!)
   - ☐ Remove pixi/rattler cache
   - ☐ Remove snakemake cache

---

## Troubleshooting

### Common Issues

**Installation fails during `pixi install`:**
- Check internet connection
- Temporarily disable antivirus
- Ensure enough disk space (~2 GB free)

**Path contains spaces error:**
- Choose a path without spaces (e.g., `C:\Users\YourName\open-tyndp` not `C:\Program Files\open-tyndp`)
- Snakemake requires paths without spaces

**Git not found error:**
- The bundled repository already includes git in the pixi environment
- Use the Start Menu shortcuts which activate the environment automatically

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
