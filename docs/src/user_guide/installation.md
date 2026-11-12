# Installation Guide

## System Requirements

* **Operating system:** Windows, macOS (Apple Silicon or Intel), or Linux.
* **CPU:** multiple cores recommended. Fastest runtimes on Apple Silicon (M-series) Macs.
* **RAM:** 16 GB minimum; 32–64 GB recommended.
* **Storage:** SSD required.

## Install the Pioneer App

Download the installer for the operating system from the [download page](https://pioneeromics.bio/download/):

| Platform | Installer |
|---|---|
| Windows | `PioneerSetup-v<version>-windows-x64.exe` |
| macOS, Apple Silicon (M1 and later) | `Pioneer-v<version>-macos-arm64.pkg` |
| macOS, Intel | `Pioneer-v<version>-macos-x64.pkg` |
| Linux | `Pioneer-v<version>-linux-x64.deb` |

Then run it:

* **Windows:** double-click the setup file and follow the wizard. It can add a desktop shortcut, and its last page offers **Launch Pioneer**. Pioneer also appears in the Start Menu.
* **macOS:** double-click the `.pkg` and follow the installer; it shows the license agreement and asks for an administrator password. The app is installed as **Pioneer v&lt;version&gt;** in `/Applications`. The installer is signed and notarized; the first launch can take about a minute while macOS verifies it.
* **Linux:** install the package with `sudo apt install ./Pioneer-v<version>-linux-x64.deb`. Pioneer appears in the applications menu.

The installer also adds the `pioneer` command for the [Command Line](@ref).

The [Quick Start Tutorial](@ref) walks through a first analysis.

## Installing Multiple Versions

Installing a different version keeps the versions already installed; reinstalling the same version replaces it. Each version gets its own app (Windows Start Menu entry, macOS app, Linux desktop entry), and the plain `pioneer` command uses the most recently installed version.

| Platform | Installed to | Run an exact version |
|---|---|---|
| Windows | `C:\Program Files\Pioneer\<version>` | use the full installation path |
| macOS | `/usr/local/lib/pioneer/<version>` | `pioneer-<version>` |
| Linux | `/opt/pioneer/<version>` | `pioneer-<version>`, or choose the default with `sudo update-alternatives --config pioneer` |

To uninstall a version on macOS, open that version's app, open **Settings**, and choose **Uninstall this version**. Settings, run history, libraries, and results are kept. On Windows, uninstall from Installed Apps; on Linux, remove the package with the package manager.

## For Developers

To run Pioneer from source in Julia (1.12 or later):

```bash
git clone https://github.com/nwamsley1/Pioneer.jl.git
cd Pioneer.jl
julia --project=.
```

```julia
julia> using Pkg; Pkg.instantiate()
julia> using Pioneer
```

The [Function Reference](@ref) lists the functions behind each step. Thermo `.raw` conversion from source requires [PioneerConverter](https://github.com/nwamsley1/PioneerConverter).
