```@meta
CurrentModule = Pioneer
```

# Pioneer

Pioneer is a fast and scalable solution for the analysis of protein mass spectrometry data. Pioneer runs through both a graphical and a command-line interface, and it scales to thousands of runs on standard desktop hardware. Pioneer includes routines for searching data-independent acquisition (DIA) experiments from Thermo Fisher, SCIEX, and Bruker instruments and for building *in silico* spectral libraries using the [Koina](https://koina.wilhelmlab.org/) interface. Prebuilt libraries for several proteomes are downloadable through Pioneer. Given a spectral library of precursor fragment ion intensities and retention time estimates, Pioneer identifies and quantifies peptides and protein groups from the library in the data.

## Key Features

- **Cross-Platform:** Pioneer runs on Linux, macOS, and Windows.
- **Multi-Vendor:** Supports Thermo Fisher, SCIEX, and Bruker instruments. See [Tested Instruments](@ref) for the instruments and mass analyzers tested with Pioneer.
- **Fast:** Pioneer searches files many times faster than data acquisition with whole-proteome *in silico* libraries, even on the latest fast-scanning TOF instruments.
- **Scalable:** Pioneer processes experiments with hundreds to thousands of runs on standard desktop hardware (32–64 GB of RAM). Search time grows linearly with the number of MS data files.
- **Deep Proteome Coverage:** Pioneer achieves high sensitivity and both quantitative precision and accuracy on benchmark datasets covering a variety of instruments and parameters.
- **Regression Testing and FDR Control:** Pioneer maintains a public suite of regression tests across diverse instruments and acquisition schemes. Each tagged release triggers searches on the benchmark datasets, which include entrapment false discovery rate and false transfer rate analyses. [Browse regression reports](https://nwamsley1.github.io/Pioneer.jl/reports/).

## Graphical User Interface

![The Pioneer app on the SearchDIA screen](assets/gui/searchdia.png)

## Get Started

- [Installation Guide](@ref): download and install Pioneer.
- [Quick Start Tutorial](@ref): a first search, step by step.
- [Parameter Configuration](@ref): what each setting does.
- [Function Reference](@ref): running Pioneer from the Julia REPL.

## Citation

Please cite the following when using Pioneer or Altimeter in published research:

> Wamsley, N. T., Wilkerson, E. M., Major, M., & Goldfarb, D. "Pioneer and Altimeter: Fast Analysis of DIA Proteomics Data Optimized for Narrow Isolation Windows." *bioRxiv* (2026). [https://www.biorxiv.org/content/10.64898/2026.02.16.706201](https://www.biorxiv.org/content/10.64898/2026.02.16.706201)

## License

©2026 Washington University

Washington University hereby grants to you a non-transferable, non-exclusive, royalty-free, non-commercial, non-clinical, not-for-use with human subjects, research license to use and copy the computer code that may be downloaded within this site (the “Software”). You agree to include this license and the above copyright notice in all copies of the Software. The Software may not be distributed, shared, or transferred to any third party. This license does not grant any rights or licenses to any other patents, copyrights, or other forms of intellectual property owned or controlled by Washington University.

If seeking a commercial license, please contact Washington University’s Office of Technology Management to request a license agreement: otm@wustl.edu.

YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED “AS IS”, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT. IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

## Help and Contact

Report problems on the [Issues page](https://github.com/nwamsley1/Pioneer.jl/issues) and ask questions or suggest features on the [Discussions page](https://github.com/nwamsley1/Pioneer.jl/discussions).

Pioneer is developed in the [Major Lab](https://majorlab.wustl.edu/) and [Goldfarb Lab](https://goldfarblab.wustl.edu/) at Washington University in St. Louis by:

- Nathan Wamsley (wamsleynathan@gmail.com)
- Dennis Goldfarb (dennis.goldfarb@wustl.edu)
