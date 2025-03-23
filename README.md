# Spectral Deconvolution via Least-Squares Fitting

This repository contains a Python script for performing spectral deconvolution using non-negative least-squares optimization. It is particularly useful in analytical chemistry, bioengineering, or spectroscopy where complex spectra need to be broken down into known pure components.

## 📌 Features

- Reads experimental and pure component spectra from CSV files
- Applies least-squares fitting with non-negativity constraints
- Outputs:
  - Optimized component coefficients for each sample
  - Reconstructed spectra
  - Goodness-of-fit metrics (SSE, RSEM)
  - Partial contributions of each component
- Ready-to-use with plotting capabilities

## 🧪 Input Files

- `YOUR_EXPERIMENTAL_dATA.csv`: Experimental spectra file with a column `wavelength` and columns for each sample.
- `PURE_SPECTRUM_REFERENCE.csv`: Pure spectra file with a column `wavelength` and columns for each reference component.

## 📤 Output Files

- `coefficients_grg.csv`: Component contribution coefficients per sample.
- `reconstructed_spectra_grg.csv`: Reconstructed spectra using the coefficients.
- `fit_statistics_grg.csv`: Goodness-of-fit statistics.
- `combined_partial_contributions_optimized.csv`: Spectral contributions from each component.

## 🚀 Usage

Just run the script with Python:

```bash
python spectral_deconvolution_final.py
```

Make sure the input files (`YOUR_EXPERIMENTAL_dATA.csv` and `PURE_SPECTRUM_REFERENCE.csv`) are in the same directory.

## 📦 Requirements

- numpy
- pandas
- scipy
- matplotlib

You can install all requirements with:

```bash
pip install -r requirements.txt
```

## 📚 License

This project is licensed under the MIT License.

## 🙌 Acknowledgements

Created with support from experimental work on spectral component analysis and deconvolution.
