# MPI Read 3D Diagnostics

MPI Read 3D is a Fortran/MPI diagnostic tool for large PadeOps 3D output fields.  It uses the same `2decomp_fft` I/O layer as PadeOps, so large fields can be read in parallel instead of loading full arrays through Python.

The active implementation is the lean direct-filename program:

```text
mpi_diag_lean.F90
```

The old keyword/budget-driven implementation has been archived in:

```text
legacy/
```

## Build

Use the unified build script:

```bash
./build_lean.sh anvil gcc
./build_lean.sh anvil impi
./build_lean.sh archer2
./build_lean.sh stampede3 gcc
./build_lean.sh stampede3 impi
```

Each local setup file in `setup/` sources the matching PadeOps environment first, then loads only the extra NetCDF module needed by this diagnostic program.  This keeps the compiler, MPI, and `2decomp_fft` stack aligned with the PadeOps build.

The build produces:

```text
build/gcc/MPIR3D_lean
build/impi/MPIR3D_lean
```

Archer2 currently has a single stack, so its default output remains
`build/MPIR3D_lean`.

You can override the default PadeOps root if needed:

```bash
PADEOPS_ROOT=/path/to/PadeOps ./build_lean.sh anvil gcc
```

## Input Model

The lean program does not use legacy budget keywords.  Fields are read by direct filenames.  Relative filenames are resolved under `path`; absolute filenames are used as-is.

Linear field composition is supported:

```text
expr {
  (1.0)(file_a.s3D)
+ (-2.0)(file_b.s3D)
+ (1.0)(file_c.s3D)
}
```

The executable reads a small namelist plus a diagnostic map:

```fortran
&SETUP
  nx       = 1330
  ny       = 466
  nz       = 700
  Lx       = 317.4603175D0
  Ly       = 55.55555555555D0
  Lz       = 55.55555555555D0
  path     = "/path/to/padeops/output"
  outdir   = "/path/to/diagnostic/output"
  driver   = "slice"
  diag_map = "examples/slice_map.diag"
/
```

Run with MPI:

```bash
mpirun -np 16 ./build/gcc/MPIR3D_lean examples/slice_input.dat
```

Use the launcher appropriate to the machine, for example `srun` on Slurm systems.

## Drivers

The current lean drivers are:

```text
ha      horizontal average
slice   2D slices and integrated slices
rms     cross-plane L2 norm profiles
profile cross-plane linear average profiles
abl     ABL height diagnostics
march   x-marched streamwise deficit reconstructions
```

For slices, all `slices {}` blocks operate on all fields:

```text
fields {
  name = u
  expr {
    (1.0)(Run05_uVel_t000900.out)
  }

  name = combo
  expr {
    (1.0)(file1.s3D)
  + (-1.0)(file2.s3D)
  }
}

slices {
  axis = y
  coords = 5.0,15.0,25.0
  integrate = 0
}
```

Wildcards in filenames are treated as a time/case dimension, not as a
summation. The program expands each wildcard pattern as a sorted list and runs
the selected driver once per case. In a composed field, each term must expand
to either one file, reused for every case, or to the same number of files as
the other wildcard terms:

```text
expr {
  (1.0)(Run05_uVel_t*.out)
+ (-1.0)(Run05_uMean.out)
}
```

Expression signs belong inside the coefficient parentheses. Use
`(-2.0)(file_b.s3D)`, not `- (2.0)(file_b.s3D)`.

This structure is intended for large batches: each field is assembled once, then every requested slice is extracted from that in-memory field.

The HA driver also supports a role-based derived field for wind speed and
wind direction:

```text
fields {
  name = wind
  derived = ws_wd

  u {
    (1.0)(Run05_uVel_t*.out)
  }

  v {
    (1.0)(Run05_vVel_t*.out)
  }
}
```

For wildcard input, the driver writes one output per match, for example
`wind_Run05_uVel_t000900_WS_HA_z.csv` and
`wind_Run05_uVel_t000900_WD_HA_z.csv`. The sorted `u` and `v` matches are
paired one-to-one, so use filename patterns that sort in timestep order.

Wind speed is the horizontal mean of pointwise `sqrt(u*u + v*v)`. Wind
direction is computed from the horizontally averaged vector,
`atan2(<v>,<u>) * 180/pi`, as a mathematical angle in degrees
counter-clockwise from +x. It is not converted to meteorological direction.

RMS-map `bounds` entries use `*` for open sides:

```text
bounds = xmin,xmax,ymin,ymax,zmin,zmax
bounds = *,*,*,*,*,*
```

The `rms` driver exports `sqrt(integral f**2 dA)`, a cross-plane L2 norm. Area normalization to form RMS is done offline.

The `rms` driver can also subtract a reference field before computing the
cross-plane L2 norm:

```text
rms {
  name = budget_combo_delta
  expr {
    (1.0)(file_a.s3D)
  + (-1.0)(file_b.s3D)
  }
  axis = x
  bounds = *,*,23.8,134.93,*,*
  subtract_reference = true
  ref_bounds = 0.0,20.0
}
```

`ref_bounds` is a two-value interval along the selected `axis`. For `axis = x`,
the reference is `f_ref(y,z)`, computed by averaging `f(x,y,z)` over the
selected `x` interval while respecting the cross-plane parts of `bounds`.
The output filename uses `_delta_rms_`, for example
`budget_combo_delta_delta_rms_x.csv`.

The `profile` driver computes linear cross-plane averages with the same
expression, axis, and bounds syntax as `rms`. For example, `axis = x` writes an
x profile of the bounded y-z average:

```text
driver = profile

profile {
  name = u_yz_avg
  expr {
    (1.0)(Run05_uVel_t000900.out)
  }
  axis = x
  bounds = *,*,10.0,45.0,0.0,55.55555556
}
```

The output is named like `u_yz_avg_avg_x.csv` and contains `x,average`.

The `march` driver reconstructs a named reference field from
`d(field)/dx = RHS`, optionally using `d(field)/dx = RHS / normalizer`. It
starts from the LES `reference` plane at `march%start`. It currently supports
`axis = x`, uses nearest grid points for `start`, `end`, and requested slice
stations, and logs the selected indices. Set `scheme = euler` or
`scheme = trapezoid`. The `reference` field is required; `normalizer` is
optional. If either named field is missing from `fields {}`, the run stops.
Every other field in `fields {}` is treated as an RHS term unless a mode lists
it under `remove`.

```text
driver = march

march {
  axis = x
  reference = deltau_les
  normalizer = uinf
  start = 200.0
  end = 220.0
  scheme = euler
  rms = true
  slices = 205.0,210.0,215.0,220.0
}

modes {
  name = G0
  remove = {
  }

  name = G1
  remove = {
    x_tiltadv_u, x_tiltadv_v, x_SGS
  }
}
```

When `rms = true`, the driver writes the same cross-plane L2 convention used by
the `rms` driver: `sqrt(integral_yz f**2 dA)`. It exports reference-field
profiles, mode reconstruction profiles, and mode error profiles, with errors
formed in 3D before the L2 profile is computed. NetCDF slice output includes
the reference field, each mode reconstruction, and each mode error.

## Examples

Copy and edit the templates in `examples/`:

```text
examples/ha_input.dat
examples/ha_map.diag
examples/ha_wind_input.dat
examples/ha_wind_map.diag
examples/slice_input.dat
examples/slice_map.diag
examples/rms_input.dat
examples/rms_map.diag
examples/profile_input.dat
examples/profile_map.diag
examples/march_map.diag
examples/abl_input.dat
examples/abl_stress_map.diag
examples/abl_inversion_map.diag
```

## Legacy

The legacy implementation remains available for reference and old workflows:

```text
legacy/mpi_diag.F90
legacy/build_Anvil.sh
legacy/build_Archer2.sh
legacy/build_Engaging.sh
legacy/build_Stampede3.sh
legacy/input.dat
legacy/run
```

Existing legacy build artifacts, if present, are also kept in `legacy/`.

New development should target `mpi_diag_lean.F90`.
