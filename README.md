# gSAM-Chem

**gSAM-Chem** is a large eddy simulation (LES) with interactive chemistry, built on the global System for Atmospheric Modeling (gSAM) version 1.8.7. It is designed to simulate atmospheric chemistry and transport at convective scales, resolving plume-scale chemistry rather than parameterizing it.

**Documentation and tutorials:** https://gsam-chem.github.io

---

## Features

- Explicit deep convection at ~66 m horizontal resolution
- Interactive gas-phase chemistry via the [Kinetic PreProcessor (KPP)](https://kpp.readthedocs.io)
- Ships with the CloudChem mechanism: 44 gas-phase + 7 photolysis reactions for tropical isoprene–NO<sub>x</sub>–OH chemistry
- Interactive surface emissions: isoprene (light- and temperature-dependent), soil NO<sub>x</sub> (precipitation-driven soil moisture), and lightning NO<sub>x</sub>
- Heterogeneous chemistry module
- Simple Land Model (SLM) available from gSAM 1.8, though the published simulations use prescribed surface latent and sensible heat fluxes from constrained variational analysis rather than SLM
- Highly parallelized (MPI); linear strong-scaling through 256 cores

---

## Prerequisites

- Fortran compiler: `ifort` (Intel Fortran Classic) or `ifx` (Intel Fortran LLVM), or `gfortran ≥ 9`
- MPI library (e.g., Intel MPI, OpenMPI, MPICH)
- NetCDF library with Fortran bindings
- `csh` shell (for the `Build` script)
- KPP 3.3.0 is only needed if you want to compile a new chemistry mechanism. The CloudChem mechanism ships with pre-generated Fortran code.

---

## Repository structure

```
gsam-chem/
├── SRC/                            # Fortran source code
│   ├── CHEMISTRY/                  # Chemistry module
│   │   ├── CLOUDCHEM/              # CloudChem KPP mechanism (.kpp, .f90)
│   │   ├── HET_CHEM/               # Heterogeneous chemistry
│   │   └── emissions.f90           # Surface and lightning emissions
│   ├── SLM/                        # Simple Land Model
│   ├── MICRO_M2005/                # Morrison two-moment microphysics
│   ├── SGS_TKE/                    # TKE-based subgrid-scale turbulence closure
│   ├── RAD_RRTM/                   # RRTM radiation
│   └── ...                         # Other dynamical core modules
├── CASES/                          # Pre-configured cases
│   └── goAMAZON_deep_convection/   # Example case (GoAmazon IOP1, wet season)
│       ├── CaseName                # Sets the case name used in output filenames
│       ├── grd                     # Vertical grid
│       ├── lsf                     # Large-scale forcing
│       ├── prm                     # Main namelist (includes &CHEMISTRY block)
│       ├── sfc                     # Surface conditions
│       └── snd                     # Initial sounding
├── SCRIPTS/                        # Post-processing and utility scripts
├── GLOBAL_DATA/                    # Shared input data (e.g., photolysis lookup tables)
├── GRIDS/                          # Grid definition files
├── UTIL/                           # Utility programs
├── Changes/                        # Changelog
├── DOC/                            # Documentation
├── OBJ/                            # Compiled object files (generated at build time)
├── Build                           # Build script (csh)
├── Makefile                        # Top-level Makefile
├── CaseName                        # Points to the active case
├── gSAM                            # Compiled executable (generated at build time)
└── gsam.run                        # Job submission script template
```

---

## Building

The build system uses a `csh` script (`Build`) that sets environment variables selecting the dynamical core modules and chemistry mechanism, then invokes the top-level `Makefile`, which sets the compiler and NetCDF paths per platform.

1. Edit `Build` to set your platform name and the desired physics/chemistry modules:

   ```csh
   setenv HOSTNAME  your-cluster-name
   setenv ADV_DIR   ADV_MPDATA
   setenv SGS_DIR   SGS_TKE
   setenv RAD_DIR   RAD_RRTM
   setenv MICRO_DIR MICRO_M2005
   setenv CHEM_MECH CLOUDCHEM
   setenv HET_CHEM  HET_CHEM
   ```

2. Edit `Makefile` to set the MPI Fortran compiler and NetCDF path for your platform (matched against `HOSTNAME`):

   ```
   FF90 = mpif90 -c -r8          # MPI Fortran compiler
   NCPATH = ${NETCDF}            # set NETCDF in your environment, or hardcode the path
   ```

3. Run the build script from the repository root:

   ```bash
   ./Build
   ```

This compiles the source in `SRC/`, placing object files in `OBJ/`, and produces the executable `gSAM` in the repository root.

**Common issue:** if the build fails with NetCDF linking errors, make sure the `NETCDF` environment variable points to your NetCDF installation prefix, or set `NCPATH` directly in the Makefile. You may also need to set `LD_LIBRARY_PATH` if NetCDF is in a non-standard location.

---

## Running a case

1. Set `CaseName` to point to your case directory (e.g., `goAMAZON_deep_convection`):

   ```bash
   echo "goAMAZON_deep_convection" > CaseName
   ```

2. Edit the namelist file (`CASES/<CaseName>/prm`) to configure the simulation. Note that the chemistry mechanism itself is chosen at **compile time** via `CHEM_MECH` in the `Build` script — the `prm` file only toggles chemistry on/off and configures submodule switches:

   ```fortran
   &PARAMETERS
     dochem = .true.,
     ...
   /

   &CHEMISTRY
     do_dry_deposition        = .true.
     do_megan_isoprene        = .true.
     do_bdsnp_no              = .true.
     do_CTG_lightning         = .true.
     do_IC_lightning          = .true.
     ...
   /
   ```

3. Submit the run (via `gsam.run` for cluster/SLURM jobs, or directly with `mpirun`):

   ```bash
   sbatch gsam.run
   # or
   mpirun -np <N> ./gSAM
   ```

Key output frequency namelist parameters: `nsave3D`, `nsave2D`, `nstat`.

Output files are written to the directory specified in `prm`, with filenames derived from `CaseName` (e.g., `goAMAZON_deep_convection_STATS.ncstat`).

---

## Adding a new chemistry mechanism

1. Generate KPP Fortran code from your `.kpp` mechanism file using KPP 3.3.0:

   ```bash
   kpp your_mechanism.kpp
   ```

2. Place the generated files in `SRC/CHEMISTRY/<YourMechanism>/`.

3. Update `SRC/CHEMISTRY/emissions.f90` as needed for mechanism-specific species.

4. Set `CHEM_MECH` to your mechanism's directory name in the `Build` script.

5. Rebuild with `./Build`.

---

## Citation

If you use gSAM-Chem in your research, please cite:

> Yoon, J.Y.S., J.A. Thornton, P.N. Blossey, M.F. Khairoutdinov, M.C. Wyant, J. Kesselmeierand, S. Kim, C.M. Nussbaumer, C. Sarkar, R. Seco, J. Shilling, J.N. Smith, A.M. Yanez-Serrano, & Turner, A.J. (submitted). gSAM–Chem: Development and Evaluation of a Large Eddy Simulation with Interactive Chemistry.

Please also cite the gSAM dynamical core:

> Khairoutdinov, M. F., Blossey, P. N., & Bretherton, C. S. (2022). Global System for Atmospheric Modeling: Model description and preliminary results. *Journal of Advances in Modeling Earth Systems*, 14(6), e2021MS002968. [https://doi.org/10.1029/2021MS002968](https://doi.org/10.1029/2021MS002968)

---

## License

See [License](https://github.com/gSAM-Chem/gSAM-Chem/blob/main/License) for details.

---

## Contact

Questions? Open an issue on [GitHub](https://github.com/gsam-chem/gsam-chem) or contact:

- **James Yoon** (lead developer): jyyoon@uw.edu
- **Alex Turner**: turneraj@uw.edu
