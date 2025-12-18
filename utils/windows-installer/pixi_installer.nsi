# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

# Alternative NSIS installer using pixi executable approach
# This installer bundles only the pixi executable and installs the environment on-demand

!include "MUI2.nsh"
!include "LogicLib.nsh"
!include "FileFunc.nsh"
!include "x64.nsh"
!include "nsDialogs.nsh"

# ------------------------------------------------------------------------------
# Installer Configuration
# ------------------------------------------------------------------------------

!define PRODUCT_NAME "Open-TYNDP"
# Version can be overridden via command line: makensis /DPRODUCT_VERSION=x.y.z
!ifndef PRODUCT_VERSION
  !define PRODUCT_VERSION "0.0.0-dev"
!endif
!define PRODUCT_PUBLISHER "Open Energy Transition"
!define PRODUCT_WEB_SITE "https://github.com/open-energy-transition/open-tyndp"
!define PRODUCT_REPO_URL "https://github.com/open-energy-transition/open-tyndp.git"
!define PIXI_ENV_NAME "open-tyndp"
!define PRODUCT_UNINST_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}"
!define PRODUCT_UNINST_ROOT_KEY "HKCU"

# Installer file name
Name "${PRODUCT_NAME} ${PRODUCT_VERSION}"
OutFile "open-tyndp-${PRODUCT_VERSION}-pixi-Windows-x86_64.exe"

# Default installation directory (user local appdata)
InstallDir "$LOCALAPPDATA\open-tyndp"

# Request user execution level (no admin required)
RequestExecutionLevel user

# Show details during installation
ShowInstDetails show
ShowUnInstDetails show

# Branding - replace "Nullsoft Install System" with OET
BrandingText "Open Energy Transition"

# Modern UI Configuration
!define MUI_ABORTWARNING
!define MUI_ICON "oet_logo.ico"
!define MUI_UNICON "oet_logo.ico"

# Branding - Open Energy Transition official colors
# Based on https://open-energy-transition.github.io/handbook/docs/Marketing/StyleGuide
!define MUI_WELCOMEFINISHPAGE_BITMAP "oet_logo.bmp"
!define MUI_HEADERIMAGE
!define MUI_HEADERIMAGE_BITMAP "${NSISDIR}\Contrib\Graphics\Header\nsis3-metro.bmp"
!define MUI_HEADERIMAGE_RIGHT

# Custom welcome page
!define MUI_WELCOMEPAGE_TITLE "Welcome to Open-TYNDP Setup"
!define MUI_WELCOMEPAGE_TITLE_3LINES
!define MUI_WELCOMEPAGE_TEXT "This wizard will guide you through the installation of Open-TYNDP, a workflow for building European energy system models.$\r$\n$\r$\nOpen-TYNDP is developed by Open Energy Transition - Accelerating the energy transition through open-source tools and data.$\r$\n$\r$\nClick Next to continue."

# Custom finish page
!define MUI_FINISHPAGE_TITLE "Installation Complete"
!define MUI_FINISHPAGE_TITLE_3LINES
!define MUI_FINISHPAGE_TEXT "Open-TYNDP has been installed successfully!$\r$\n$\r$\nUse the Start Menu shortcuts under 'Open-TYNDP' to launch PowerShell or Command Prompt with the environment activated.$\r$\n$\r$\nFor documentation and support, visit:"
!define MUI_FINISHPAGE_LINK "https://github.com/open-energy-transition/open-tyndp"
!define MUI_FINISHPAGE_LINK_LOCATION "https://github.com/open-energy-transition/open-tyndp"

# ------------------------------------------------------------------------------
# Macros
# ------------------------------------------------------------------------------

# Locate Windows binaries (cmd.exe, powershell.exe)
# Based on https://github.com/conda/constructor/blob/main/constructor/nsis/main.nsi.tmpl
!macro FindWindowsBinaries
    # Find cmd.exe and powershell.exe
    ReadEnvStr $R0 SystemRoot
    ReadEnvStr $R1 windir
    ${If} ${FileExists} "$R0"
        StrCpy $CMD_EXE "$R0\System32\cmd.exe"
        StrCpy $POWERSHELL_EXE "$R0\System32\WindowsPowerShell\v1.0\powershell.exe"
    ${ElseIf} ${FileExists} "$R1"
        StrCpy $CMD_EXE "$R1\System32\cmd.exe"
        StrCpy $POWERSHELL_EXE "$R1\System32\WindowsPowerShell\v1.0\powershell.exe"
    ${Else}
        # Cross our fingers binaries are in PATH
        StrCpy $CMD_EXE "cmd.exe"
        StrCpy $POWERSHELL_EXE "powershell.exe"
    ${EndIf}
!macroend

# ------------------------------------------------------------------------------
# Variables
# ------------------------------------------------------------------------------

Var /GLOBAL REPO_DIR
Var /GLOBAL CLONE_REPO
Var /GLOBAL PIXI_EXE
Var /GLOBAL CMD_EXE
Var /GLOBAL POWERSHELL_EXE

# Custom page variables
Var /GLOBAL Dialog
Var /GLOBAL Label
Var /GLOBAL DirRequest
Var /GLOBAL DirBrowse
Var /GLOBAL CheckBox
Var /GLOBAL CheckBoxState

# Uninstaller page variables
Var /GLOBAL UnDialog
Var /GLOBAL UnCheckBoxRepo
Var /GLOBAL UnCheckBoxPixiCache
Var /GLOBAL UnCheckBoxSnakemakeCache
Var /GLOBAL UnRemoveRepo
Var /GLOBAL UnRemovePixiCache
Var /GLOBAL UnRemoveSnakemakeCache

# ------------------------------------------------------------------------------
# Installer Pages
# ------------------------------------------------------------------------------

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "..\..\LICENSES\MIT.txt"

# Custom page for repository setup
Page custom RepositoryPageCreate RepositoryPageLeave

!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

# Uninstaller pages
UninstPage custom un.UninstallOptionsPageCreate un.UninstallOptionsPageLeave
!insertmacro MUI_UNPAGE_INSTFILES

# Language
!insertmacro MUI_LANGUAGE "English"

# ------------------------------------------------------------------------------
# Functions for Custom Repository Page
# ------------------------------------------------------------------------------

Function RepositoryPageCreate
    !insertmacro MUI_HEADER_TEXT "Repository Setup" "Choose where to clone and install Open-TYNDP"

    nsDialogs::Create 1018
    Pop $Dialog

    ${If} $Dialog == error
        Abort
    ${EndIf}

    # Info text
    ${NSD_CreateLabel} 0 0 100% 24u "The Open-TYNDP repository contains the workflow scripts and configuration files.$\r$\n$\r$\nThis installer will clone the repository and use pixi to install the environment automatically.$\r$\n$\r$\nNote: Internet connection required for downloading packages from conda-forge."
    Pop $Label

    # Checkbox to enable/disable cloning
    ${NSD_CreateCheckbox} 0 28u 100% 12u "Clone repository and install environment automatically"
    Pop $CheckBox
    ${NSD_SetState} $CheckBox ${BST_CHECKED}
    ${NSD_OnClick} $CheckBox RepositoryCheckBoxClick

    # Directory label
    ${NSD_CreateLabel} 0 48u 100% 12u "Repository directory:"
    Pop $Label

    # Directory text field
    ${NSD_CreateDirRequest} 0 62u 85% 12u "$PROFILE\open-tyndp"
    Pop $DirRequest

    # Browse button
    ${NSD_CreateBrowseButton} 86% 62u 14% 12u "Browse..."
    Pop $DirBrowse
    ${NSD_OnClick} $DirBrowse RepositoryBrowseClick

    # Installation location info
    ${NSD_CreateLabel} 0 82u 100% 24u "The pixi executable will be installed to:$\r$\n%LOCALAPPDATA%\open-tyndp\pixi.exe$\r$\n$\r$\nThe environment will be installed inside the repository directory using pixi."
    Pop $Label

    nsDialogs::Show
FunctionEnd

Function RepositoryCheckBoxClick
    Pop $0
    ${NSD_GetState} $0 $CheckBoxState
    ${If} $CheckBoxState == ${BST_CHECKED}
        EnableWindow $DirRequest 1
        EnableWindow $DirBrowse 1
    ${Else}
        EnableWindow $DirRequest 0
        EnableWindow $DirBrowse 0
    ${EndIf}
FunctionEnd

Function RepositoryBrowseClick
    Pop $0
    nsDialogs::SelectFolderDialog "Select Repository Directory" "$PROFILE"
    Pop $0
    ${If} $0 != error
        ${NSD_SetText} $DirRequest "$0\open-tyndp"
    ${EndIf}
FunctionEnd

Function RepositoryPageLeave
    # Get the checkbox state
    ${NSD_GetState} $CheckBox $CLONE_REPO

    # Get the directory
    ${NSD_GetText} $DirRequest $REPO_DIR

    # Store in registry for later use by launcher scripts
    WriteRegStr HKCU "Software\${PRODUCT_NAME}" "RepositoryPath" "$REPO_DIR"

    ${If} $CLONE_REPO == ${BST_CHECKED}
        # Check if directory already exists
        ${If} ${FileExists} "$REPO_DIR\.git"
            MessageBox MB_YESNO|MB_ICONQUESTION \
                "A Git repository already exists at:$\n$REPO_DIR$\n$\nSkip setup and use existing repository?" \
                /SD IDYES \
                IDYES skip_clone
            # User chose to not skip, abort installation
            Abort
        ${EndIf}

        ${If} ${FileExists} "$REPO_DIR\*.*"
            MessageBox MB_YESNO|MB_ICONQUESTION \
                "The directory already exists and is not empty:$\n$REPO_DIR$\n$\nDo you want to continue anyway?" \
                /SD IDNO \
                IDYES continue_anyway
            Abort
            continue_anyway:
        ${EndIf}

        skip_clone:
    ${EndIf}
FunctionEnd

# ------------------------------------------------------------------------------
# Installer Section
# ------------------------------------------------------------------------------

Section "Install" SecInstall
    SetOutPath "$INSTDIR"

    # Locate Windows binaries
    !insertmacro FindWindowsBinaries

    DetailPrint "Installing ${PRODUCT_NAME} ${PRODUCT_VERSION}..."

    # Create the installation directory if it doesn't exist
    CreateDirectory "$INSTDIR"

    # Copy the pixi executable
    # NOTE: You need to download pixi.exe for Windows first
    DetailPrint "Installing pixi executable..."
    File "pixi.exe"

    # Set pixi path
    StrCpy $PIXI_EXE "$INSTDIR\pixi.exe"

    # Verify pixi works
    DetailPrint "Verifying pixi installation..."
    nsExec::ExecToLog '"$CMD_EXE" /D /C "$\"$PIXI_EXE$\" --version"'
    Pop $0
    ${If} $0 != 0
        MessageBox MB_OK|MB_ICONEXCLAMATION "Failed to verify pixi executable. Error code: $0"
        Abort "Installation failed"
    ${EndIf}

    # Clone and setup repository if requested
    ${If} $CLONE_REPO == ${BST_CHECKED}
        ${If} ${FileExists} "$REPO_DIR\.git"
            DetailPrint "Repository already exists at $REPO_DIR"

            # Check if it needs environment installation
            ${If} ${FileExists} "$REPO_DIR\.pixi\envs\default"
                DetailPrint "Environment already installed, skipping..."
                Goto skip_env_install
            ${Else}
                DetailPrint "Environment not found, will install..."
                Goto do_env_install
            ${EndIf}
        ${Else}
            DetailPrint "Cloning repository to $REPO_DIR..."
            DetailPrint "This may take several minutes depending on your internet connection..."

            # Create parent directory if needed
            CreateDirectory "$REPO_DIR"

            # Clone using pixi exec git
            DetailPrint "Cloning repository using pixi exec git..."
            # Note: --verbose provides text output without progress bar control characters
            nsExec::ExecToLog '"$CMD_EXE" /D /C "$\"$PIXI_EXE$\" exec git clone --verbose $\"${PRODUCT_REPO_URL}$\" $\"$REPO_DIR$\""'
            Pop $0
            ${If} $0 != 0
                MessageBox MB_YESNO|MB_ICONEXCLAMATION \
                    "Failed to clone repository (error code: $0).$\n$\nDo you want to continue installation anyway?" \
                    /SD IDYES \
                    IDYES continue_without_clone
                Abort "Installation cancelled"
                continue_without_clone:
                DetailPrint "Repository clone failed, but continuing installation..."
                Goto skip_env_install
            ${Else}
                DetailPrint "Repository cloned successfully!"
            ${EndIf}
        ${EndIf}

        do_env_install:
        # Install the environment using pixi
        DetailPrint "Installing environment using pixi..."
        DetailPrint "This will download ~500-800 MB of packages from conda-forge..."
        DetailPrint "This may take 5-15 minutes depending on your internet connection..."

        # Change to repository directory and run pixi install
        # Use -v for info-level logging to show download/install progress
        SetOutPath "$REPO_DIR"
        nsExec::ExecToLog '"$CMD_EXE" /D /C "$\"$PIXI_EXE$\" install -e ${PIXI_ENV_NAME} -v"'
        Pop $0
        ${If} $0 != 0
            MessageBox MB_YESNO|MB_ICONEXCLAMATION \
                "Failed to install environment (error code: $0).$\n$\nDo you want to continue installation anyway?" \
                /SD IDYES \
                IDYES continue_without_env
            Abort "Installation cancelled"
            continue_without_env:
            DetailPrint "Environment installation failed, but continuing..."
        ${Else}
            DetailPrint "Environment installed successfully!"
        ${EndIf}

        skip_env_install:
    ${EndIf}

    # Create Start Menu shortcuts directly (no .bat files needed)
    DetailPrint "Creating Start Menu shortcuts..."
    CreateDirectory "$SMPROGRAMS\${PRODUCT_NAME}"

    # PowerShell shortcut - launches PowerShell and executes pixi shell inside it
    CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME}\Open-TYNDP PowerShell.lnk" \
        "$POWERSHELL_EXE" \
        "-ExecutionPolicy ByPass -Command $\"Set-Location '$REPO_DIR'; & '$INSTDIR\pixi.exe' shell -e ${PIXI_ENV_NAME}$\"" \
        "$POWERSHELL_EXE" \
        0 \
        SW_SHOWNORMAL \
        "" \
        "Launch PowerShell with Open-TYNDP environment"

    # CMD shortcut - launches cmd and executes pixi shell inside it
    CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME}\Open-TYNDP Command Prompt.lnk" \
        "$CMD_EXE" \
        "/C $\"cd /d $\"$REPO_DIR$\" && $\"$INSTDIR\pixi.exe$\" shell -e ${PIXI_ENV_NAME}$\"" \
        "$CMD_EXE" \
        0 \
        SW_SHOWNORMAL \
        "" \
        "Launch Command Prompt with Open-TYNDP environment"

    # Uninstaller shortcut
    CreateShortCut "$SMPROGRAMS\${PRODUCT_NAME}\Uninstall Open-TYNDP.lnk" \
        "$INSTDIR\Uninstall.exe" \
        "" \
        "$INSTDIR\Uninstall.exe" \
        0 \
        SW_SHOWNORMAL \
        "" \
        "Uninstall Open-TYNDP"

    # Write the uninstaller
    DetailPrint "Creating uninstaller..."
    WriteUninstaller "$INSTDIR\Uninstall.exe"

    # Write registry keys for Add/Remove Programs
    WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayName" "${PRODUCT_NAME}"
    WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayVersion" "${PRODUCT_VERSION}"
    WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "Publisher" "${PRODUCT_PUBLISHER}"
    WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "URLInfoAbout" "${PRODUCT_WEB_SITE}"
    WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "UninstallString" "$INSTDIR\Uninstall.exe"
    WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "InstallLocation" "$INSTDIR"
    WriteRegDWORD ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "NoModify" 1
    WriteRegDWORD ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "NoRepair" 1

    # Write estimated size (pixi.exe is small, ~20 MB)
    WriteRegDWORD ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "EstimatedSize" 20480

    DetailPrint "Installation completed successfully!"

SectionEnd

# ------------------------------------------------------------------------------
# Uninstaller Page Functions
# ------------------------------------------------------------------------------

Function un.UninstallOptionsPageCreate
    !insertmacro MUI_HEADER_TEXT "Uninstall Options" "Choose what to remove"

    nsDialogs::Create 1018
    Pop $UnDialog

    ${If} $UnDialog == error
        Abort
    ${EndIf}

    # Info text
    ${NSD_CreateLabel} 0 0 100% 24u "Select which components you want to remove:$\r$\n$\r$\nWARNING: Removing the repository will delete all your work and configuration files!"
    Pop $Label

    # Read repository path from registry
    ReadRegStr $REPO_DIR HKCU "Software\${PRODUCT_NAME}" "RepositoryPath"

    # Repository checkbox
    ${NSD_CreateCheckbox} 0 30u 100% 12u "Remove repository directory and all work files"
    Pop $UnCheckBoxRepo
    ${If} $REPO_DIR != ""
    ${AndIf} ${FileExists} "$REPO_DIR"
        ${NSD_SetState} $UnCheckBoxRepo ${BST_UNCHECKED}
    ${Else}
        EnableWindow $UnCheckBoxRepo 0
    ${EndIf}

    # Repository path label
    ${If} $REPO_DIR != ""
    ${AndIf} ${FileExists} "$REPO_DIR"
        ${NSD_CreateLabel} 20u 44u 95% 12u "Location: $REPO_DIR"
        Pop $Label
    ${EndIf}

    # Pixi cache checkbox
    ${NSD_CreateCheckbox} 0 62u 100% 12u "Remove pixi/rattler cache directory"
    Pop $UnCheckBoxPixiCache
    ${NSD_SetState} $UnCheckBoxPixiCache ${BST_UNCHECKED}

    # Pixi cache path label
    ${NSD_CreateLabel} 20u 76u 95% 12u "Location: %LOCALAPPDATA%\rattler\cache"
    Pop $Label

    # Snakemake cache checkbox
    ${NSD_CreateCheckbox} 0 94u 100% 12u "Remove snakemake-pypsa-eur cache directory"
    Pop $UnCheckBoxSnakemakeCache
    ${NSD_SetState} $UnCheckBoxSnakemakeCache ${BST_UNCHECKED}

    # Snakemake cache path label
    ${NSD_CreateLabel} 20u 108u 95% 12u "Location: %LOCALAPPDATA%\snakemake-pypsa-eur"
    Pop $Label

    nsDialogs::Show
FunctionEnd

Function un.UninstallOptionsPageLeave
    # Get checkbox states
    ${NSD_GetState} $UnCheckBoxRepo $UnRemoveRepo
    ${NSD_GetState} $UnCheckBoxPixiCache $UnRemovePixiCache
    ${NSD_GetState} $UnCheckBoxSnakemakeCache $UnRemoveSnakemakeCache
FunctionEnd

# ------------------------------------------------------------------------------
# Uninstaller Section
# ------------------------------------------------------------------------------

Section "Uninstall"
    # Locate Windows binaries
    !insertmacro FindWindowsBinaries

    DetailPrint "Removing ${PRODUCT_NAME}..."

    # Remove repository if requested
    ${If} $UnRemoveRepo == ${BST_CHECKED}
        ReadRegStr $REPO_DIR HKCU "Software\${PRODUCT_NAME}" "RepositoryPath"
        ${If} $REPO_DIR != ""
            # Remove .pixi environment directory first
            Push "$REPO_DIR\.pixi"
            Push "pixi environment directory"
            Call un.RemoveDirectory

            # Remove repository directory
            Push "$REPO_DIR"
            Push "repository directory"
            Call un.RemoveDirectory
        ${EndIf}
    ${EndIf}

    # Remove pixi cache if requested
    ${If} $UnRemovePixiCache == ${BST_CHECKED}
        Push "$LOCALAPPDATA\rattler\cache"
        Push "pixi/rattler cache"
        Call un.RemoveDirectory
    ${EndIf}

    # Remove snakemake cache if requested
    ${If} $UnRemoveSnakemakeCache == ${BST_CHECKED}
        Push "$LOCALAPPDATA\snakemake-pypsa-eur"
        Push "snakemake-pypsa-eur cache"
        Call un.RemoveDirectory
    ${EndIf}

    # Remove Start Menu shortcuts
    DetailPrint "Removing Start Menu shortcuts..."
    Delete "$SMPROGRAMS\${PRODUCT_NAME}\Open-TYNDP PowerShell.lnk"
    Delete "$SMPROGRAMS\${PRODUCT_NAME}\Open-TYNDP Command Prompt.lnk"
    Delete "$SMPROGRAMS\${PRODUCT_NAME}\Uninstall Open-TYNDP.lnk"
    RMDir "$SMPROGRAMS\${PRODUCT_NAME}"

    # Remove installation directory
    DetailPrint "Removing installation files..."
    Delete "$INSTDIR\pixi.exe"
    Delete "$INSTDIR\Uninstall.exe"

    # Remove any remaining files and the installation directory
    RMDir "$INSTDIR"

    # Remove registry keys
    DeleteRegKey ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}"
    DeleteRegKey HKCU "Software\${PRODUCT_NAME}"

    DetailPrint "Uninstallation completed!"

SectionEnd

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

# Function to remove a directory quickly using cmd RMDIR followed by NSIS RMDir
# Parameters:
#   $0 - Directory path to remove
#   $1 - Description of what's being removed
Function un.RemoveDirectory
    Pop $1  # Description
    Pop $0  # Directory path

    ${If} ${FileExists} "$0"
        DetailPrint "Removing $1: $0"
        # Use cmd RMDIR first (faster for large directories)
        nsExec::Exec '"$CMD_EXE" /D /C RMDIR /Q /S "$0"'
        # Follow up with NSIS RMDir to ensure complete removal
        RMDir /r "$0"
    ${Else}
        DetailPrint "$1 not found at: $0"
    ${EndIf}
FunctionEnd

Function .onInit
    # Set default repository directory
    StrCpy $REPO_DIR "$PROFILE\open-tyndp"
    StrCpy $CLONE_REPO ${BST_CHECKED}

    # Check if already installed
    ReadRegStr $R0 ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "UninstallString"
    ${If} $R0 != ""
        MessageBox MB_YESNO|MB_ICONQUESTION \
            "${PRODUCT_NAME} is already installed.$\n$\nDo you want to uninstall the existing version first?" \
            /SD IDYES \
            IDYES uninst
        Abort

        uninst:
            # Run the uninstaller
            ExecWait '$R0 _?=$INSTDIR'
            # Check if uninstallation was successful
            ${If} ${FileExists} "$INSTDIR\Uninstall.exe"
                MessageBox MB_OK|MB_ICONEXCLAMATION "Uninstallation failed or was cancelled. Please uninstall manually before installing again."
                Abort
            ${EndIf}
    ${EndIf}
FunctionEnd
