# Installation Guide

## System Requirements
* **Julia**: 1.10 or higher
* **CPU**: Multiple cores and threads recommended. Increasing the number of threads reduces computation time.
* **RAM**: >=16GB recommended. RAM availability should exceed the spectral library size by at least 4GB. For searching against the yeast proteome, as little as 6-8 GB may suffice.  
* **Storage**: SSD recommended. Available disk space at least double the total size of the `.arrow` formmated raw files to search. The `.arrow` files are usually ~1/2 the size of the vendor files. 
* **Operating System**: Windows, Linux, or macOS

## Installation Steps
### 1. Install Julia

1. Download and install Julia 1.10 or higher from [julia.org](https://julialang.org/downloads/) 
2. Verify installation by opening the Julia REPL and checking the version:
```@julia
julia> versioninfo()
Julia Version 1.10.0
...
```

### 2. Install Pioneer.jl
#### Standard Instructions 
1. Start the Julia REPL and install 
```@julia
julia> using Pioneer
[ Info: Precompiling Pioneer [585db54e-1982-41c9-ae07-9e9eb56c7d61]
┌ Info: Skipping precompilation due to precompilable error. Importing Pioneer [585db54e-1982-41c9-ae07-9e9eb56c7d61].
└ exception = Error when precompiling module, potentially caused by a __precompile__(false) declaration in the module.
```
2. Test your installation. From the REPL, enter package mode (press ]), then:
```@julia
(@v1.11) pkg> test Pioneer
```

#### Developeer Instructions
To install Pioneer for development, clone the GitHub repository and work in the
dedicated development environment. This keeps production dependencies clean
while providing tools like Documenter, PProf, Revise, and PackageCompiler.
1. Open Command Prompt (Windows) or Terminal (MacOSX/Linux) and clone the repository:
   ```@julia
   git clone https://github.com/nwamsley1/Pioneer.jl.git
   cd Pioneer.jl
   ```
2. Start Julia with the development environment and enter package mode (press `]`):
   ```@julia
   julia --project=dev
   pkg> develop ./
   ```
3. Return from package mode (backspace) and verify installation:
   ```@julia
   julia> using Pioneer
   ```
4. Test your installation. Enter package mode (press `]`), then:
   ```@julia
   pkg> test Pioneer
   ```

!!! note "'note'"
    On the first attempt ```using Pioneer``` requires an internet connection, and it may take several minutes to download and install dependencies.

!!! note "'note'"
    The best way to develop Pioneer is with the help of the [`Revise` package](https://timholy.github.io/Revise.jl/stable/).

### 3. Install PioneerConverter

Pioneer requires Thermo RAW files to be converted to Apache Arrow format. The cross-platform [PioneerConverter](https://github.com/nwamsley1/PioneerConverter/releases/tag/v0.1.0) package does the conversion. 

1. Download .NET 8.0 SDK and .NET Runtime 8.0 [here](https://dotnet.microsoft.com/en-us/download/dotnet/8.0)

2. Download and install PioneerConverter
   * Using Release (recommended): 
      * Download and decompress the latest release of [PioneerConverter](https://github.com/nwamsley1/PioneerConverter/releases/tag/v0.1.0)

   * Build locally: 
      * Clone the repository
         ```
         git clone https://github.com/nwamsley1/PioneerConverter.git
         ```
      * Navitage inside of the PioneerConverter directory and run: ```dotnet build -c Release```
3. Test your installation:
   * Using the Command Prompt (windows) or Terminal (MaxOS/Linus), navigate inside of the PioneerConverter directory.
      * Windows: `cmd bin\Release\net8.0\PioneerConverter.exe \path\to\test\raw\file.raw -b 10000`
      * Linux/macOS: `bin/Release/net8.0/PioneerConverter /path/to/test/raw/file.raw -b 10000`

!!! note "'note'"
    Pioneer.jl has experimental support for converting .mzML formatted SCIEX TOF data to the .arrow format via the `convertMzML` method. 

### Next Steps

After installation:
1. Follow the [Quick Start Tutorial](@ref)
2. Review [Parameter Configuration](@ref "Parameter Configuration")