# Computational certificate for common-zero contours

This repository contains the computational artifact accompanying the paper:

Agasthya Chidambaram, *Divisible formal weight enumerators: Chebyshev determinants, extremal bounds, and zeros of Duursma zeta polynomials*.

The program verifies the finite algebraic inequalities used in the computer-assisted part of Theorem C.

## Contents

- `verify_common_zero_contour.py` — verification program
- `certified_output.txt` — complete output from the archived verification run
- `requirements.txt` — pinned Python-package requirements
- `environment.txt` — software environment used for the archived run
- `SHA256SUMS` — SHA-256 checksums of the computational files
- `CITATION.cff` — software citation metadata
- `LICENSE` — software license

## Requirements

The exact environment used for the archived verification is recorded in
`environment.txt`. The required Python packages are pinned in
`requirements.txt`.

Python 3.10 or later is required. The verifier must be run without Python's
`-O` optimization option because the certificate conditions are enforced by
assertions.

The certificate inequalities use exact `fractions.Fraction` arithmetic.
SymPy is used separately to count real algebraic roots in the rational
isolating intervals specified in the paper.

## Reproduction on Windows

Open PowerShell in the repository directory and create a clean virtual
environment:

```powershell
py -m venv .venv
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r .\requirements.txt
```

Run the verification:

```powershell
python .\verify_common_zero_contour.py
```

Do not run the verifier with `python -O`, since optimized execution disables
assertions.

A successful execution terminates with exit code zero and ends with:

```text
All certificate checks passed.
```

The complete expected output is stored in `certified_output.txt`.

## Comparing the output

To save a new run and compare it with the archived output, use:

```powershell
$NewOutput = & python .\verify_common_zero_contour.py 2>&1

$NewOutput |
    ForEach-Object { $_.ToString() } |
    Set-Content -Encoding utf8 .\new_output.txt

$Difference = Compare-Object `
    (Get-Content .\certified_output.txt) `
    (Get-Content .\new_output.txt)

if ($Difference) {
    $Difference
    Write-Host "The outputs differ."
}
else {
    Write-Host "The outputs agree exactly."
}
```

The temporary file may then be removed:

```powershell
Remove-Item .\new_output.txt
```

## Integrity verification

The SHA-256 checksums of the archived computational files are recorded in
`SHA256SUMS`.

To recompute them in PowerShell, run:

```powershell
$FilesToHash = @(
    ".\verify_common_zero_contour.py",
    ".\certified_output.txt",
    ".\requirements.txt",
    ".\environment.txt"
)

Get-FileHash $FilesToHash -Algorithm SHA256 |
    Sort-Object Path |
    ForEach-Object {
        "{0}  {1}" -f $_.Hash.ToLower(), (Split-Path $_.Path -Leaf)
    }
```

The resulting values should agree exactly with `SHA256SUMS`.

## Scope of the computation

The script performs the following finite checks:

1. SymPy root counts confirm the stated rational isolating intervals.
2. Exact rational interval arithmetic verifies the required positivity and
   upper-bound inequalities.
3. Bernstein-basis coefficients certify the required polynomial positivity.
4. The ring-specific base estimates are checked for each of the cases
   \(r=4,5,6\).

The mathematical justification for these checks and their role in the proof
are given in the accompanying paper.

## Repository

The development repository is available at:

https://github.com/cookie579/verify-common-zero-contour

The permanent archived release and DOI are provided by Zenodo. The
version-specific Zenodo release should be used when citing the software.

## Citation

Citation metadata are supplied in `CITATION.cff`. Please cite the
version-specific Zenodo record associated with the release used in your work.

## Author

Agasthya Chidambaram  
Faculty of Mathematics, University of Waterloo  
Waterloo, Ontario, Canada  
a3chidam@uwaterloo.ca

## License

The source code is released under the MIT License. See `LICENSE` for the
complete license terms.