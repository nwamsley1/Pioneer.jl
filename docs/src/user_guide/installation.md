# Installation Guide

## System Requirements
* **CPU**: Multiple cores and threads recommended. Increasing the number of threads reduces computation time.
* **RAM**: >=16GB recommended. RAM availability should exceed the spectral library size by at least 4GB. For searching against the yeast proteome, as little as 6-8 GB may suffice.
* **Storage**: SSD recommended. Available disk space at least double the total size of the `.arrow` formatted raw files to search. The `.arrow` files are usually ~1/2 the size of the vendor files.
* **Operating System**: Windows, Linux, or macOS

## Installation

### End-User Installation
1. Download the installer for your operating system from the [releases page](https://github.com/nwamsley1/Pioneer.jl/releases):
   * **Windows** – `PioneerSetup-*.exe` installer
   * **macOS** – `.pkg` installer (signed and notarized); separate builds for Intel and Apple Silicon
   * **Linux** – `.deb` package
2. Run the installer. It places a `pioneer` executable on your `PATH`.
3. On first launch:
   * **macOS** – Gatekeeper verifies the binary and the first run can take about a minute. Zipped binaries require manual Gatekeeper approval and are not recommended.
4. Verify the installation:
   ```bash
   pioneer --help
   ```
5. Use `--threads N` to control the number of worker threads:
   ```bash
   pioneer search --threads 8 ...
   ```

### Installing Multiple Versions

Installer releases are versioned independently. Reinstalling the same Pioneer
version replaces that version, while installing another version keeps both:

| Platform | Versioned installation | Unqualified `pioneer` command |
|---|---|---|
| Windows | `C:\Program Files\Pioneer\<version>` | Most recently installed version; takes effect in new terminals |
| macOS | `/usr/local/lib/pioneer/<version>` | Most recently installed version via `/usr/local/lib/pioneer/current` |
| Linux | `/opt/pioneer/<version>` | Highest installed release selected by `update-alternatives` |

Windows lists each version separately in Installed Apps and creates a
version-labelled Start Menu shortcut. macOS installs version-labelled apps in
`/Applications`. Linux installs version-labelled desktop entries and also
provides commands such as `pioneer-2.1.0` for selecting an exact version.

The Windows setup offers an optional, version-labelled desktop shortcut. Its
success page also provides **Launch Pioneer**, which starts the newly installed
GUI. For quiet deployment, pass `CreateDesktopShortcut=1` to the setup
executable to request the shortcut.

On Linux, inspect or override the active version with:

```bash
sudo update-alternatives --config pioneer
```

On macOS, exact versions are available as `/usr/local/bin/pioneer-<version>`.
On Windows, invoke an exact version using its full installation path.

To uninstall a macOS version, open that version-labelled Pioneer app, open
Settings, and choose **Uninstall this version**. macOS asks for an administrator
password, then removes that app and its matching command-line tools. Pioneer
settings, run history, spectral libraries, and analysis results are left alone.
If another version remains installed, the unqualified `pioneer` command is
automatically redirected to the most recently installed remaining version.

### Docker
Run Pioneer in a container without installing dependencies.

1. Pull the prebuilt image:
   ```bash
   docker pull dennisgoldfarb/pioneer:latest
   ```
2. Execute Pioneer inside the container, mounting a host directory (e.g. the current directory) to access data:
   ```bash
   docker run --rm -it -v $(pwd):/work dennisgoldfarb/pioneer:latest pioneer --help
   ```
   Replace `pioneer --help` with any subcommand.
3. To build the image locally using the included `Dockerfile`:
   ```bash
   docker build -t pioneer .
   ```

### Development Setup
To work on Pioneer itself, set up a local development environment.

1. Install **Julia 1.10** or higher from [julialang.org](https://julialang.org/downloads/) (only required for development; compiled releases bundle their own runtime).
2. Clone the repository:
   ```bash
   git clone https://github.com/nwamsley1/Pioneer.jl.git
   cd Pioneer.jl
   ```
3. Start Julia in the development environment and activate the project:
   ```julia
   julia --project=dev
   pkg> develop ./
   ```
4. In the Julia REPL load Revise and Pioneer:
   ```julia
   julia> using Revise, Pioneer
   ```
5. Install [PioneerConverter](https://github.com/nwamsley1/PioneerConverter) to convert Thermo RAW files to Arrow format.
6. Call the main functions directly, e.g.
   ```julia
   # Option 1: Single FASTA directory (backward compatible)
   params = GetBuildLibParams(out_dir, lib_name, fasta_dir)
   BuildSpecLib(params)
   
   # Option 2: Flexible input - files and/or directories
   params = GetBuildLibParams(out_dir, lib_name, 
       ["/path/to/dir1", "/path/to/file.fasta", "/path/to/dir2"])
   BuildSpecLib(params)
   params = GetSearchParams("library.poin", "ms_data", "results")
   SearchDIA(params)
   ```

| Subcommand       | Julia function   |
|------------------|------------------|
| `params-predict` | `GetBuildLibParams` |
| `predict`        | `BuildSpecLib`     |
| `params-search`  | `GetSearchParams`  |
| `search`         | `SearchDIA`        |
| `convert-raw`    | `PioneerConverter` |
| `convert-mzml`   | `convertMzML`      |

!!! note
    `Revise` enables hot reloading of code during development.

### PioneerConverter
Detailed installation and usage instructions for PioneerConverter are available in its [documentation](https://github.com/nwamsley1/PioneerConverter).
Common `convert-raw` options are `--output-dir`, `--skip-existing`, `--concurrent-files`, and `--threads-per-file`. For full option details, run `pioneer convert-raw --help`.

### Next Steps

After installation:
1. Follow the [Quick Start Tutorial](@ref).
2. Generate parameter files with `pioneer params-predict` or `pioneer params-search`,
   then edit them according to [Parameter Configuration](@ref "Parameter Configuration").
