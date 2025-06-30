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

| Option        | Description                                                                       |
| ------------- | --------------------------------------------------------------------------------- |
| `-e`          | Error threshold (default: `1e-10`)                                                |
| `-p`          | Working precision in decimal digits                                               |
| `-P`          | Number of decimal digits to write to file (in HDF5 files, written as `double`)    |
| `-t`          | Number of threads                                                                 |
| `-a`          | **Abscissa** in Durbin’s method, used as a safety margin when `-e` is set         |
| `-w`          | Weight for distributing error between consistency and truncation (default: `0.5`) |
| `-N`          | Maximum number of Laplace summands in Durbin's method (Laplace inversion)         |
| `-n`          | Number of Hankel summands in Dini’s expansion (Hankel inversion)                  |
| `-kahan`      | Enable **Kahan summation**                                                        |
| `-ignorecstr` | Ignore **CSTRs** in error estimation                                              |
| `-model`      | Input model file (`.h5` or `.xml`)                                                |
| `-o`          | Output file (`.csv` or `.h5`). Optional for HDF5; defaults to input file          |

---

## Notes on Error Estimation

In some configurations (e.g., **cyclic systems**), CASEMA **does not support error estimation**.

* Specifying the error threshold (`-e`) in such cases will result in an exception.
* Instead, the user must specify:

  * The **abscissa** using `-a`
  * The **number of Hankel summands** using `-n`

---

Let me know if you'd like this reformatted into a PDF, HTML page, or auto-generated documentation site (e.g., using MkDocs or Docusaurus).
