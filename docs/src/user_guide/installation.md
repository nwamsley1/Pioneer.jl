# Installation Guide

## System Requirements
* **Julia**: 1.10 or higher
* **CPU**: Multiple cores and threads recommended. Increasing the number of threads reduces computation time.
* **RAM**: >=16GB recommended. RAM availability should exceed the spectral library size by at least 4GB. For searching against the yeast proteome, as little as 6-8 GB may suffice.  
* **Storage**: SSD recommended. Available disk space at least double the total size of the `.arrow` formmated raw files to search. The `.arrow` files are usually ~1/2 the size of the vendor files. 
* **Operating System**: Windows, Linux, or macOS

## Installation

### End-User Installation
1. Download the installer for your operating system from the [releases page](https://github.com/nwamsley1/Pioneer.jl/releases).
2. Run the installer. It places a `pioneer` executable on your `PATH`.
3. On first launch:
   * **macOS** – Gatekeeper verifies the binary and the first run can take about a minute. Zipped binaries require manual Gatekeeper approval and are not recommended.
4. Verify the installation:
   ```bash
   pioneer --help
   ```

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

1. Install Julia 1.10 or higher from [julia.org](https://julialang.org/downloads/).
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

!!! note "'note'"
    `Revise` enables hot reloading of code during development.

### PioneerConverter
Detailed installation and usage instructions for PioneerConverter are available in its [documentation](https://github.com/nwamsley1/PioneerConverter).

### Next Steps

After installation:
1. Follow the [Quick Start Tutorial](@ref).
2. Generate parameter files with `pioneer params-predict` or `pioneer params-search`,
   then edit them according to [Parameter Configuration](@ref "Parameter Configuration").
