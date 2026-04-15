---

# CADET-Semi-Analytic (CASEMA) Interface Documentation

The interface of **CADET-Semi-Analytic (CASEMA)** follows the **CADET file format v50000**, as documented at [https://cadet.github.io/](https://cadet.github.io/). This format is implemented, for example, in **CADET-Core version 5.0.0**.

However, there are key differences between CASEMA and CADET-Core due to differences in simulatable models and numerical methods.

---

## Key Differences from CADET-Core

### 1. Supported Models

* CASEMA **only supports single-component systems**.
* The only available binding model is the **linear binding model**.

### 2. 2D GRM Outlet Handling

* **Outlet ports in 2D GRM models** represent **volume-averaged values**.
* To compare solutions, a recommended approach is to define an outlet for each radial zone, assigning all finite volume (FV) cells of each zone to their respective outlet.

### 3. Output Behavior

* Only **chromatograms (outlet solutions)** are always written for all unit operations.
* The following output fields **do not exist** and are **implicitly disabled**:

  * `write_solution_inlet`
  * `write_solution_bulk`
  * `write_solution_particle`
  * `write_solution_solid`
* The field `write_solution_outlet` is **implicitly true** for all units.

### 4. Parameter Sensitivities

* Parameter sensitivities are **deprecated**, but still **available in version 1**.

---

## Numerical Settings and Configuration

CASEMA does **not** use the `solver` or `discretization` groups from the CADET input file.
Instead, the following program options can be provided to control numerical methods, tolerances, and I/O.

### Program Options

- `-model <file>` (model input file required)
  - Path to the input model file. Accepted file formats are `.h5` (HDF5) and `.xml`, both have to follow the CADET input format.
- `-m <file>` (model input file required)
  - Path to the input model file. Accepted file formats are `.h5` (HDF5) and `.xml`, both must follow the CADET input format.

- `-o <file>` (solution output file optional)
- `-o <file>` (solution output file, optional)
  - Path to the output file. Supported formats are `.csv` and `.h5`.
  - For `.h5`: program options are saved as meta data. Loss of arbitrary precision output: Numeric values are written as doubles. Working precision is not affected.
  - For `.csv`: arbitrary precision solution, the number formatting uses the `-P` option;
  - Optional for `.h5` inputs (in that case CASEMA will default the output to the input filename). For `.xml` inputs specify an output path explicitly if you want a `.h5` or `.csv` result.
  - HDF5 (`.h5`) behaviour: numeric solution arrays are written as IEEE `double` (64-bit). CASEMA writes meta data, including the program options into the file.
  - CSV (`.csv`) behaviour: arbitrary precision solution, number formatting is based on the `-P` option;
  - Optional for `.h5` inputs (in that case CASEMA will default the output to the input filename). For `.xml` inputs specify an output path explicitly if you want a `.h5` or `.csv` result.

- `-e <threshold>` (error threshold, default: `1e-10`)
  - Relative error threshold used by CASEMA's adaptive inversion/error-estimation framework. The solver guarantees that the numerical error (consistency + Laplace truncation) does not exceed this value.
  - A smaller `-e` requests a more accurate solution. To meet a smaller threshold CASEMA will generally increase the number of Laplace summands and may increase internal precision, which increases runtime and memory usage.
  - See "Notes on Error Estimation" for exclusions and further information.

- `-w <weight>` (consistency / truncation weight, default: `0.5`)
  - A weight in `[0,1]` that distributes the allowable error between consistency error (numerical discretization/rounding) and truncation error (series/inversion truncation). `0` assigns all budget to truncation error, `1` to consistency error.
  - Adjusts how CASEMA balances the two error sources when attempting to meet `-e`. Tweaking `-w` can change the number of summands used vs. the internal solver settings; it affects numerical trade-offs but not the physical model.
  - Optional, recommended for advanced users only. The default `0.5` is suitable for most problems.

- `-a <abscissa>` (abscissa for Durbin's Laplace inversion)
  - The abscissa is the real shift of the Bromwich contour used by Durbin's Laplace inversion. It acts as a safety margin controlling the location of the inversion path in the complex plane.
  - Choosing a suitable `-a` affects convergence of the Laplace inversion. If `-a` is too small the inversion may converge slowly or be unstable; if `-a` is unnecessarily large truncation errors can increase. `-a` does not change the underlying physical model, only the numerical inversion behavior.
  - Required when error estimation is not available (e.g., cyclic systems). When `-e` is set and error estimation is active, CASEMA computes an internal `-a` if none is provided; providing `-a` manually overrides that choice and is an advanced tuning option.

- `-N <max-summands>` (maximum Durbin Laplace summands)
  - Upper limit for the number of Laplace summands that Durbin's method may use during inversion.
  - A higher `-N` allows the inversion to use more terms to reach the requested accuracy; if too small the method may stop early and the result will be dominated by truncation error. Increasing `-N` increases CPU time and, for very large values, memory usage.
  - Increase `-N` only if CASEMA reports that the default cap is reached and accuracy is insufficient. If you set `-e` and the algorithm requires more summands than `-N` allows, a warning or error will be emitted.

- `-n <hankel-summands>` (Hankel/Dini summands)
  - Required only for 2D models (2DGRM), otherwise ignored
  - Number of Hankel summands used in Dini's expansion for Hankel inversion steps.
  - Directly controls truncation error in Hankel expansions: larger `-n` reduces truncation error and produces a more accurate radial transform at the cost of more computation.

- `-p <digits>` (working precision in decimal digits)
  - Sets the arithmetic precision used internally (number of decimal digits) when an arbitrary-precision backend is available (MPFR). If unset, CASEMA uses the current MPFR default precision (visible in the CLI as the default `-p`).
  - Effect on solution: Higher working precision reduces round-off error in ill-conditioned transforms (many summands, very small `-e`) and can be necessary to attain extreme tolerances. Increasing `-p` can dramatically increase runtime and memory usage.
  - Only set when instability is observed or if its impossible to reach `-e` with double precision.
  - Optional, CASEMA uses standard IEEE double precision by default.

- `-P <digits>` (output precision)
  - Number of decimal digits written to text output (CSV). If `-P` is zero/unspecified the program uses the working precision `-p` as the default for text formatting. `-P` does not affect binary HDF5 storage (HDF5 arrays are doubles).

- `-t <threads>` (number of threads)
  - Maximum number of worker threads CASEMA should use for parallel parts of the computation.

- `-ignorecstr` (ignore CSTRs in error estimation)
  - Specifies whether CSTRs are omitted from the automatic error-estimation procedure.
  - Influences how the error budget is distributed across the system; ignoring CSTRs may loosen the estimated error bound and reduce the required number of summands. It does not change the computed deterministic solution other than via the different error-driven numerical settings.
  - Optional, recommended for advanced users. This option is useful if you know CSTRs in your model are fast-mixing and should not dominate the error estimate. Otherwise leave it disabled to get conservative error estimates.

- `-kahan` (enable Kahan summation)
  - Enables Kahan compensated summation for long floating-point sums.
  - Reduces rounding bias in long accumulations (e.g., when summing many small terms), which can improve numerical accuracy without changing the algorithm. There is a small runtime cost.
  - Optional, recommended when using many summands (large `-N`/`-n`) or when accumulation errors are observed.

---

## Notes on Error Estimation

CASEMA supports two mutually exclusive modes for controlling Laplace/Hankel inversion accuracy:

- Automatic (error-driven): provide `-e <threshold>`, optionally together with a weight `-w` (0 < `w` < 1).
CASEMA guarantees error bounds and computes a safe abscissa and the per-unit number of Laplace summands required to meet the requested relative error.
If you additionally specify `-N`, it acts as an upper cap for the automatically chosen summands.
An exception is raised when guaranteed error bounds are not available.

- Manual (abscissa + summands): provide `-a <abscissa>` and `-N <max-summands>`.
This mode is **required** when guaranteed error bounds are not available.
Guaranteed error bounds are not available in:
- cyclic systems
- if any unit operation has initial conditions $> 0$
- CSTRs with variable volume
- GRM with surface diffusion and kinetic adsorption
- 2D models (2DGRM). Here, the user must additionally supply a positive `-n <hankel-summands>`.
