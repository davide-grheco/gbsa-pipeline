---
title: Configuration Reference
---

# Configuration Reference

All pipeline settings are declared in a single TOML file. Every section is
optional except `[system]`.

GROMACS MDP parameter names use **hyphens** in `.mdp` files but **underscores**
in the TOML config (e.g., `ref-t` → `ref_t`).

---

## `[system]`

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `protein` | path | **yes** | Path to the protein PDB file |
| `ligand` | path | no | Path to the ligand SDF file (3-D conformer required) |
| `extra_ff_files` | list[path] | no | Extra OpenMM ForceField XML files (e.g., metal parameters) |
| `net_charge` | int | no | Formal charge of the ligand (auto-detected when omitted) |

---

## `[forcefield]`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `protein_ff` | str | `"ff14SB"` | Protein force field. Options: `"ff14SB"`, `"ff19SB"`, `"ff99SB"` |
| `ligand_ff` | str | `"gaff2"` | Ligand force field. Options: `"gaff"`, `"gaff2"` |
| `charge_method` | str | `"am1bcc"` | Partial charge method. Options: `"am1bcc"`, `"nagl"`, `"espaloma-am1bcc"` |

---

## `[solvation]`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `water_model` | str | `"tip3p"` | Water model. Options: `"tip3p"`, `"tip4p"`, `"tip5p"`, `"spc"`, `"spce"` |
| `box_shape` | str | `"truncated_octahedron"` | Box shape. Options: `"truncated_octahedron"`, `"cubic"` |
| `padding` | float | `null` | Distance (nm) from solute to box edge. Takes precedence over `box_size` |
| `box_size` | float | `8.0` | Absolute box edge length (nm). Used when `padding` is not set |
| `ion_concentration` | float | `0.15` | Salt concentration (mol/L) |
| `neutralize` | bool | `true` | Add counter-ions to neutralize the system |

Either `padding` or `box_size` must be provided.

---

## `[minimization]`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `nsteps` | int | `10000` | Maximum number of minimization steps |
| `emtol` | float | `10.0` | Convergence criterion: max force (kJ mol⁻¹ nm⁻¹) |

---

## `[equilibration]`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `simulation_time_ps` | float | `500.0` | NVT heating duration (ps), ramping from 0 K to 300 K |

---

## `[md]` — Production MD (GROMACS MDP parameters)

The `[md]` section accepts any field of
[`GromacsParams`][gbsa_pipeline.change_defaults.GromacsParams].
Field names use underscores; they map to hyphenated GROMACS MDP keys.

### Integrator

| TOML field | MDP key | Type | Default | Description |
|------------|---------|------|---------|-------------|
| `integrator` | `integrator` | str | `"md"` | Integration algorithm: `"md"` (leapfrog), `"md-vv"`, `"sd"`, `"bd"` |
| `dt` | `dt` | float | `0.001` | Time step (ps) |
| `nsteps` | `nsteps` | int | `500` | Number of MD steps |
| `tinit` | `tinit` | float | `0.0` | Start time (ps) |

### Output control

| TOML field | MDP key | Type | Default | Description |
|------------|---------|------|---------|-------------|
| `nstlog` | `nstlog` | int | `500` | Steps between log entries |
| `nstenergy` | `nstenergy` | int | `500` | Steps between energy writes |
| `nstxout_compressed` | `nstxout-compressed` | int | `500` | Steps between compressed trajectory frames |

### Thermostat

| TOML field | MDP key | Type | Default | Description |
|------------|---------|------|---------|-------------|
| `tcoupl` | `tcoupl` | str | `"no"` | Thermostat: `"no"`, `"berendsen"`, `"nose-hoover"`, `"v-rescale"`, `"andersen"` |
| `tc_grps` | `tc-grps` | str | `"System"` | Temperature coupling group(s) |
| `tau_t` | `tau-t` | float | `0.1` | Temperature coupling time constant (ps) |
| `ref_t` | `ref-t` | float | `300.0` | Reference temperature (K) |
| `nhchainlength` | `nhchainlength` | int | `10` | Nosé-Hoover chain length |

### Barostat

| TOML field | MDP key | Type | Default | Description |
|------------|---------|------|---------|-------------|
| `pcoupl` | `pcoupl` | str | `"no"` | Barostat: `"no"`, `"Berendsen"`, `"Parrinello-Rahman"`, `"C-rescale"`, `"MTTK"` |
| `pcoupltype` | `pcoupltype` | str | `"isotropic"` | Coupling geometry: `"isotropic"`, `"semiisotropic"`, `"anisotropic"`, `"surface-tension"` |
| `tau_p` | `tau-p` | float | `2.0` | Pressure coupling time constant (ps) |
| `ref_p` | `ref-p` | float | `1.0` | Reference pressure (bar) |
| `compressibility` | `compressibility` | float | `4.5e-5` | Isothermal compressibility (bar⁻¹) |

### Electrostatics & VdW

| TOML field | MDP key | Type | Default | Description |
|------------|---------|------|---------|-------------|
| `coulombtype` | `coulombtype` | str | `"PME"` | Electrostatics method |
| `rcoulomb` | `rcoulomb` | float | `1.2` | Coulomb cutoff (nm) |
| `vdw_type` | `vdw-type` | str | `"Cut-off"` | VdW interaction type |
| `rvdw` | `rvdw` | float | `1.2` | VdW cutoff (nm) |
| `rvdw_switch` | `rvdw-switch` | float | `1.0` | VdW switch start (nm) |

### Constraints

| TOML field | MDP key | Type | Default | Description |
|------------|---------|------|---------|-------------|
| `constraints` | `constraints` | str | `"none"` | Constraint type: `"none"`, `"h-bonds"`, `"all-bonds"`, `"h-angles"`, `"all-angles"` |
| `constraint_algorithm` | `constraint-algorithm` | str | `"LINCS"` | Solver: `"LINCS"` or `"SHAKE"` |
| `lincs_order` | `lincs-order` | int | `4` | LINCS expansion order |

### Velocity generation

| TOML field | MDP key | Type | Default | Description |
|------------|---------|------|---------|-------------|
| `gen_vel` | `gen-vel` | str | `"no"` | Generate initial velocities: `"yes"` or `"no"` |
| `gen_temp` | `gen-temp` | float | `300.0` | Temperature for velocity generation (K) |
| `gen_seed` | `gen-seed` | int | `-1` | Random seed (`-1` = use system clock) |
