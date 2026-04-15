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

CASEMA does **not** use the `solver` or `discretization` groups from the CADET input file. Instead, numerical options must be provided **via command-line arguments**.

### Command-Line Options

The following command-line options control numerical methods, tolerances and I/O. Each option is described with its effect on the computed solution, sensible defaults, and guidance when it must or should be set (or when it will be ignored).

- `-model <file>` (required)
  - Description: Path to the input model file. Accepted formats are `.h5` (HDF5) and `.xml` (CADET XML).
  - Effect: The model file fully defines the simulated system. CASEMA cannot run without a valid model file.
  - When to set: Always; this is required.

- `-o <file>` (optional)
  - Description: Path to the output file. Supported formats are `.csv` and `.h5`.
  - Effect: Controls where computed chromatograms and any auxiliary output are written. For `.csv` the number formatting uses the `-P` option; for `.h5` numeric values are written as doubles and metadata is preserved where possible.
  - When to set: Optional for `.h5` inputs (in that case CASEMA will default the output to the input filename). For `.xml` inputs specify an output path explicitly if you want a `.h5` or `.csv` result.

- `-e <threshold>` (error threshold, default: `1e-10`)
  - Description: Relative error threshold used by CASEMA's adaptive inversion/error-estimation framework. The solver tries to guarantee that the numerical error (consistency + truncation) does not exceed this value.
  - Effect on solution: A smaller `-e` requests a more accurate solution. To meet a smaller threshold CASEMA will generally increase the number of inversion summands and internal precision, which increases runtime and memory usage. If `-e` is too large results may be less accurate.
  - When to set: Recommended when you need a reproducible accuracy target (for example, comparing with reference data). Do not set in configurations where CASEMA does not support error estimation (see "Notes on Error Estimation" below) — in those cases providing `-e` will raise an exception and you must use `-a` and `-n` instead.

- `-a <abscissa>` (abscissa for Durbin's Laplace inversion)
  - Description: The abscissa is the real shift of the Bromwich contour used by Durbin's Laplace inversion. It acts as a safety margin controlling the location of the inversion path in the complex plane.
  - Effect on solution: Choosing a suitable `-a` affects convergence of the Laplace inversion. If `-a` is too small the inversion may converge slowly or be unstable; if `-a` is unnecessarily large truncation errors can increase. `-a` does not change the underlying physical model, only the numerical inversion behavior.
  - When to set: Required when error estimation is not available (e.g., cyclic systems). When `-e` is set and error estimation is active, CASEMA computes an internal `-a` if none is provided; providing `-a` manually overrides that choice and is an advanced tuning option.

- `-N <max-summands>` (maximum Durbin Laplace summands)
  - Description: Upper limit for the number of Laplace summands that Durbin's method may use during inversion.
  - Effect on solution: A higher `-N` allows the inversion to use more terms to reach the requested accuracy; if too small the method may stop early and the result will be dominated by truncation error. Increasing `-N` increases CPU time and, for very large values, memory usage.
  - When to set: Increase `-N` only if CASEMA reports that the default cap is reached and accuracy is insufficient. If you set `-e` and the algorithm requires more summands than `-N` allows, a warning or error will be emitted.

- `-n <hankel-summands>` (Hankel/Dini summands)
  - Description: Number of Hankel summands used in Dini's expansion for Hankel inversion steps (used in radial / cylinder transforms and some GRM kernels).
  - Effect on solution: `-n` directly controls truncation error in Hankel expansions: larger `-n` reduces truncation error and produces a more accurate radial transform at the cost of more computation. If `-n` is too small the radial approximation will be visibly inaccurate in the solution.
  - When to set: Required to be explicitly provided in configurations where automatic error estimation is disabled or unsupported (see "Notes on Error Estimation"). Otherwise CASEMA will choose a safe default; explicitly increasing `-n` is recommended when validating convergence for radial problems.

- `-p <digits>` (working precision in decimal digits)
  - Description: Sets the arithmetic precision used internally (number of decimal digits) when an arbitrary-precision backend is available. If not specified, CASEMA uses standard IEEE double precision by default.
  - Effect on solution: Higher working precision reduces round-off error in ill-conditioned transforms (many summands, very small `-e`) and can be necessary to attain extreme tolerances. Increasing `-p` can dramatically increase runtime and memory usage.
  - When to set: Only set when you observe instability or inability to reach `-e` with double precision, or when the model contains very stiff or ill-conditioned transforms. For typical problems double precision is sufficient and `-p` can be left unset.

- `-P <digits>` (output precision / formatting)
  - Description: Number of decimal digits written to text output (CSV). For HDF5 output numeric values are written as doubles; `-P` only controls how many digits are printed when exporting to human-readable formats.
  - Effect on solution: Does not affect computation — only controls formatting of output values.
  - When to set: Set when you need a specific decimal output resolution for downstream processing or visual inspection. Leave unset to use a sensible default formatting.

- `-t <threads>` (number of threads)
  - Description: Maximum number of worker threads CASEMA should use for parallel parts of the computation.
  - Effect on solution: Primarily affects performance (wall-clock time). In principle multi-threading can change the order of floating-point operations and therefore produce small numerical differences due to rounding; such differences are normally negligible. If your environment or build uses a single-threaded numerical backend this option will be ignored.
  - When to set: Set to the number of CPU cores you want CASEMA to utilize. Leave unset to let CASEMA choose a sensible default (often the machine's core count).

- `-w <weight>` (consistency / truncation weight, default: `0.5`)
  - Description: A number in `[0,1]` that distributes the allowable error between consistency error (numerical discretization/rounding) and truncation error (series/inversion truncation). `0` assigns all budget to truncation error, `1` to consistency error.
  - Effect on solution: Adjusts how CASEMA balances the two error sources when attempting to meet `-e`. Tweaking `-w` can change the number of summands used vs. the internal solver settings; it affects numerical trade-offs but not the physical model.
  - When to set: Advanced users only. The default `0.5` is suitable for most problems.

- `-kahan` (enable Kahan summation)
  - Description: Enables Kahan compensated summation for long floating-point sums.
  - Effect on solution: Reduces rounding bias in long accumulations (e.g., when summing many small terms), which can improve numerical accuracy without changing the algorithm. There is a small runtime cost.
  - When to set: Recommended when using many summands (large `-N`/`-n`) or when you observe small systematic errors that look like accumulation error. If your platform lacks a compensated summation implementation the option will be ignored.

- `-ignorecstr` (ignore CSTRs in error estimation)
  - Description: When enabled, continuous stirred-tank reactors (CSTRs) are omitted from the automatic error-estimation procedure.
  - Effect on solution: This influences how the error budget is distributed across the system; ignoring CSTRs may loosen the estimated error bound and reduce the required number of summands. It does not change the computed deterministic solution other than via the different error-driven numerical settings.
  - When to set: Use with care. This option is useful if you know CSTRs in your model are fast-mixing and should not dominate the error estimate. Otherwise leave it disabled to get conservative error estimates.


---

## Notes on Error Estimation

In some configurations (e.g., **cyclic systems**), CASEMA **does not support error estimation**.

* Specifying the error threshold (`-e`) in such cases will result in an exception.
* Instead, the user must specify:

  * The **abscissa** using `-a`
  * The **number of Hankel summands** using `-n`

