# gSAM-Chem

**gSAM-Chem** is a large eddy simulation (LES) with interactive chemistry, built on the global System for Atmospheric Modeling ([gSAM](http://rossby.msrc.sunysb.edu/~marat/SAM.html)) version 1.8. It is designed to simulate atmospheric chemistry and transport at convective scales, resolving plume-scale chemistry rather than parameterizing it.

**Documentation and tutorials:** https://gsam-chem.github.io

---

## Features

- Explicit deep convection at ~50 m horizontal resolution
- Interactive gas-phase chemistry via the [Kinetic PreProcessor (KPP)](https://kpp.readthedocs.io)
- Ships with the **CloudChem mechanism**: 44 gas-phase + 7 photolysis reactions for tropical isoprene–NO<sub>x</sub>–OH chemistry
- Interactive surface emissions: isoprene (light- and temperature-dependent), soil NO<sub>x</sub> (precipitation-driven soil moisture), and lightning NO<sub>x</sub>
- Heterogeneous chemistry module
- Simple Land Model (SLM) from gSAM 1.8
- Highly parallelized (MPI); superlinear strong-scaling through 128+ cores

---

## Prerequisites

- Fortran compiler: `ifort` (Intel Fortran Classic) or `ifx` (Intel Fortran LLVM)
- MPI library (e.g., Intel MPI, OpenMPI)
- NetCDF library with Fortran bindings
- `csh` shell (for the `Build` script)
- KPP is only needed if you want to compile a new chemistry mechanism. The CloudChem mechanism ships with pre-generated Fortran code.

---

## Repository structure

```
gsam-chem/
├── SRC/                        # Fortran source code
│   ├── CHEMISTRY/              # Chemistry module
│   │   ├── CLOUDCHEM/          # CloudChem KPP mechanism (.kpp, .f90)
│   │   ├── HET_CHEM/           # Heterogeneous chemistry
│   │   └── emissions.f90       # Surface and lightning emissions
│   ├── SLM/                    # Simple Land Model
│   ├── MICRO_M2005/            # Morrison two-moment microphysics
│   ├── SGS_TKE/                # TKE-based subgrid-scale turbulence closure
│   ├── RAD_RRTM/               # RRTM radiation
│   └── ...                     # Other dynamical core modules
├── CASES/
│   └── goAMAZON_deep_convection/   # Example case (GoAmazon IOP2)
│       ├── CaseName            # Sets the case name used in output filenames
│       ├── grd                 # Vertical grid
│       ├── lsf                 # Large-scale forcing
│       ├── prm                 # Main namelist (includes &CHEMISTRY block)
│       ├── sfc                 # Surface conditions
│       └── snd                 # Initial sounding
├── SCRIPTS/                    # Post-processing and utility scripts
├── GLOBAL_DATA/                # Shared input data (e.g., photolysis lookup tables)
├── GRIDS/                      # Grid definition files
├── UTIL/                       # Utility programs
├── DOC/                        # Documentation
├── OBJ/                        # Compiled object files (generated at build time)
├── Build                       # Build script (csh)
├── Makefile                    # Top-level Makefile
├── CaseName                    # Points to the active case
├── gSAM                        # Compiled executable (generated at build time)
└── gsam.run                    # Job submission script template
```

---

## Building

The build system uses a `csh` script that sets environment variables and invokes the top-level `Makefile`.

1. Edit `Build` to set your compiler, MPI paths, and NetCDF paths:

   ```csh
   setenv FC ifort          # or ifx
   setenv MPIFC mpif90
   setenv NETCDF /path/to/netcdf
   ```

2. Run the build script from the repository root:

   ```bash
   ./Build
   ```

This compiles the source in `SRC/`, placing object files in `OBJ/`, and produces the executable `gSAM` in the repository root.

---

## Running a case

1. Set `CaseName` to point to your case directory (e.g., `CASES/goAMAZON_deep_convection`):

   ```bash
   echo "goAMAZON_deep_convection" > CaseName
   ```

2. Edit the namelist file (`CASES/<CaseName>/prm`) to configure the simulation. Chemistry options are set in the `&CHEMISTRY` block within `prm`:

   ```fortran
   &CHEMISTRY
     dochem        = .true.
     mech_name     = 'CLOUDCHEM'
     ...
   /
   ```

3. Run the model (or submit via `gsam.run` for cluster jobs):

   ```bash
   mpirun -np <N> ./gSAM
   ```

Key output frequency namelist parameters: `nsave3D`, `nsave2D`, `nstat`.

Output files are written to the directory specified in `prm`, with filenames derived from `CaseName` (e.g., `goAMAZON_deep_convection_STATS.ncstat`).

---

## Adding a new chemistry mechanism

1. Generate KPP Fortran code from your `.kpp` mechanism file:

   ```bash
   kpp your_mechanism.kpp
   ```

2. Place the generated files in `SRC/CHEMISTRY/<YourMechanism>/`.

3. Update `SRC/CHEMISTRY/emissions.f90` as needed for mechanism-specific species.

4. Rebuild with `./Build`.

---

## Citation

If you use gSAM-Chem in your research, please cite:

> Yoon, J. Y. S., Wyant, M., Blossey, P., Thornton, J. A., & Turner, A. J. (in preparation). *gSAM-Chem: Implementing chemistry into a large eddy simulation.*

Please also cite the gSAM dynamical core:

> Khairoutdinov, M. F. & Randall, D. A. (2003). Cloud resolving modeling of the ARM summer 1997 IOP: Model formulation, results, uncertainties, and sensitivities. *Journal of the Atmospheric Sciences*, 60(4), 607–625. https://doi.org/10.1175/1520-0469(2003)060<0607:CRMOTA>2.0.CO;2

---

## License

The gSAM-Chem chemistry module and case files are released under the MIT License. See [LICENSE](LICENSE) for details. The gSAM dynamical core is subject to its own license terms — see the gSAM documentation.

---

## Contact

Questions? Open an issue on [GitHub](https://github.com/gsam-chem/gsam-chem) or contact:

- **James Yoon** (lead developer): jyyoon@uw.edu
- **Alex Turner**: turneraj@uw.edu
