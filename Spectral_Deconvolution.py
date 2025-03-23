"""
Spectral Deconvolution via Least-Squares Fitting

Author: [Sunil Tiwari]
Date: [March-23-2025]
License: MIT
Repository: https://github.com/yourusername/spectral-deconvolution

Description:
------------
This script performs spectral deconvolution using a least-squares optimization approach, nnls model. It is specifically tailored for spectroscopic data analysis, where 
experimental spectra are expressed as linear combinations of reference spectra (often from pure compounds).

The script:
1. Loads experimental and standard (pure compound) spectra.
2. Performs non-negative least-squares fitting for each sample.
3. Calculates reconstructed spectra and goodness-of-fit metrics.
4. Outputs deconvoluted coefficients and reconstruction results.
5. (Optional) Generates plots to visualize the fitting quality.

Input:
------
- Your experimrntal.csv: Experimental spectra with wavelength in the first column.
- Standard.csv: Pure reference spectra with wavelength in the first column.

Output:
-------
- coefficients_grg.csv: Optimized contribution of each component to each sample.
- reconstructed_spectra_grg.csv: Fitted spectra reconstructed from coefficients.
- fit_statistics_grg.csv: Statistical metrics including SSE and RSEM for each sample.
- combined_partial_contributions_optimized.csv: Spectral contributions of each component.
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# 1) LEAST-SQUARES FIT #
def run_least_squares_fitting():
    """
    Performs the steps from 'this script':
    1. Loads 'Your experimental csv file' & 'pure spectra csv file'.
    2. Minimize SSE subject to non-negativity for each sample.
    3. Save:
       - coefficients_grg.csv ##Coeffieients it gets from deconvoluted spectra
       - reconstructed_spectra_grg.csv  #Sum of individual components in the spectra added 
       - fit_statistics_grg.csv  ##SSE and RSEM of fit
       - combined_partial_contributions_optimized.csv (needed for subsequent steps)
    4. Also returns the wavelength array for optional plotting.
    """
    # -- Load Data --
    samples_df = pd.read_csv("SLs_evolved_strans_chl_again.csv") ##Input your experimental file name
    pure_df    = pd.read_csv("pure_spectra_final.csv")  ##Input your pure component file

    wavelengths = samples_df['wavelength'].values	##Make sure your first column heading is wavelength not Wavelength
    D = samples_df.drop(columns=['wavelength']).values   # shape: (n_wavelengths, n_samples)
    ST = pure_df.drop(columns=['wavelength']).values     # shape: (n_wavelengths, n_components)

    n_samples = D.shape[1]
    n_components = ST.shape[1]
    n_wavelengths = D.shape[0]
    print(f"Loaded {n_samples} samples, {n_components} components, {n_wavelengths} wavelengths.")

    # -- Objective Function (SSE) --
    def objective_function(coefficients, ST, measured_spectrum):
        reconstructed_spectrum = ST @ coefficients
        residuals = measured_spectrum - reconstructed_spectrum
        return np.sum(residuals**2)
    # -- Non-Negativity Constraints --
    def create_constraints(n_components):
        constraints = []
        # Each coefficient >= 0
        for i in range(n_components):
            constraints.append({'type': 'ineq', 'fun': lambda x, i=i: x[i]})
        return constraints

    # -- Solve for Each Sample --
    coefficients_results = []
    reconstructed_spectra = []
    sse_list = []
    rmse_list = []

    for i in range(n_samples):
        measured_spectrum = D[:, i]
        x0 = np.full(n_components, 0.1)  # initial guess
        constraints = create_constraints(n_components)
        result = minimize(
            objective_function,
            x0,
            args=(ST, measured_spectrum),
            method='SLSQP',
            constraints=constraints,
            bounds=[(0, None)] * n_components
        )

        if result.success:
            print(f"Sample {i+1}: Optimization succeeded.")
        else:
            print(f"Sample {i+1}: Optimization failed. Reason: {result.message}")

        # Store solutions
        coefficients_results.append(result.x)
        rec_spectrum = ST @ result.x
        reconstructed_spectra.append(rec_spectrum)

        # SSE / RMSE
        residuals = measured_spectrum - rec_spectrum
        sse = np.sum(residuals**2)
        sse_list.append(sse)
        rmse = np.sqrt(sse / len(measured_spectrum))
        rmse_list.append(rmse)

    # -- Convert to arrays --
    coefficients_results = np.array(coefficients_results)   # shape: (n_samples, n_components)
    reconstructed_spectra = np.array(reconstructed_spectra) # shape: (n_samples, n_wavelengths)
    # -- Save Results --
    # 1. Coefficients to CSV

    coeff_cols = pure_df.columns[1:]  # name of the pure components
    sample_names = samples_df.columns[1:]
    coefficients_df = pd.DataFrame(coefficients_results, columns=coeff_cols, index=sample_names)
    coefficients_df.to_csv("coefficients_grg.csv")
    print("Saved coefficients to 'coefficients_grg.csv'.")

    # 2. Reconstructed Spectra to CSV
    #    -> shape must be (n_wavelengths, n_samples) for a typical DataFrame format
    rec_spectra_df = pd.DataFrame(reconstructed_spectra.T, columns=sample_names)
    rec_spectra_df.insert(0, 'wavelength', wavelengths)
    rec_spectra_df.to_csv("reconstructed_spectra_grg.csv", index=False)
    print("Saved reconstructed spectra to 'reconstructed_spectra_grg.csv'.")

    # 3. Fit statistics (SSE, RMSE) to CSV
    stats_df = pd.DataFrame({
        'Sample': sample_names,
        'SSE': sse_list,
        'RMSE': rmse_list
    })
    stats_df.to_csv("fit_statistics_grg.csv", index=False)
    print("Saved fit statistics to 'fit_statistics_grg.csv'.")

    # 4. Combined partial contributions
    #    - For each sample, scale each component's spectrum by its coefficient
    partial_contributions_list = []
    base_df = pd.DataFrame({'wavelength': wavelengths})
    for i, sample_name in enumerate(sample_names):
        # scale columns in ST by each coefficient for sample i
        scaled = ST * coefficients_results[i]  # shape: (n_wavelengths, n_components)
        # rename columns to reflect sample + component
        col_names = [f"{sample_name}_{comp}" for comp in coeff_cols]
        temp_df = pd.DataFrame(scaled, columns=col_names)
        partial_contributions_list.append(temp_df)

    # combine side by side
    all_partial_contributions_df = pd.concat([base_df] + partial_contributions_list, axis=1)
    all_partial_contributions_df.to_csv("combined_partial_contributions_optimized.csv", index=False)
    print("Saved partial contributions to 'combined_partial_contributions_optimized.csv'.")

    # Optional quick plot for Sample 1
    plt.figure(figsize=(10, 6))
    plt.plot(wavelengths, D[:, 0], label="Measured Spectrum (Sample 1)", linewidth=2)
    plt.plot(wavelengths, reconstructed_spectra[0], label="Reconstructed Spectrum (Sample 1)", linestyle="--")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Absorbance")
    plt.legend()
    plt.title("Measured vs Reconstructed Spectrum (Sample 1)")
    plt.show()
    return wavelengths  # or any other data you'd like to reuse for further analysis

###############################
# 2) EXTRACT PEAK ABSORBANCES #
###############################
def extract_peak_absorbances(

    input_csv="combined_partial_contributions_optimized.csv",
    output_csv="extracted_peak_absorbances.csv"
):
    """
    Reads a CSV of deconvoluted spectra (columns like 'Sample_Component'),
    extracts the absorbance at specific wavelengths for each known pigment
    (parsed from the last underscore in the column name), and writes them out.
    """
    # 1. Dictionary of pigment -> target wavelength (Here you can change the pigments lambda max based on pigments needed or change for respective pigments)
    peak_wavelengths = {
        "B,B-carotene": 450,
        "Myxoxanthophyll": 478,
        "Zeaxanthin": 450,
        "Echinenone": 458,
        "Chlorophyll-a": 665
    }

    df = pd.read_csv(input_csv)
    if "wavelength" not in df.columns:
        raise ValueError("Input CSV must contain a 'wavelength' column.")
    results_list = []
    for col in df.columns:
        if col.lower() == "wavelength":
            continue
        # e.g. "Sample1_B,B-carotene"
        parts = col.split("_")
        pigment_name = parts[-1]  # last token
        sample_id = "_".join(parts[:-1])
        if pigment_name in peak_wavelengths:
            target_wl = peak_wavelengths[pigment_name]
            row_match = df.loc[df["wavelength"] == target_wl]

            if len(row_match) == 1:
                absorbance_value = row_match[col].values[0]
            else:
                absorbance_value = None
            results_list.append({
                "Original_Column": col,
                "Sample_ID": sample_id,
                "Pigment": pigment_name,
                "Target_Wavelength": target_wl,
                "Extracted_Absorbance": absorbance_value
            })
    out_df = pd.DataFrame(results_list)
    out_df.to_csv(output_csv, index=False)
    print(f"Extraction complete. Results saved to '{output_csv}'.")

################################
# 3) CONCENTRATION CALCULATION #
################################
def compute_concentration_ug_per_mL(
    input_csv="extracted_peak_absorbances.csv",
    output_csv="with_concentrations.csv"
):
    """
    Reads the extracted peak absorbances CSV, applies Beer–Lambert law for 1 cm pathlength
    to compute concentration in µg/mL for each pigment, and saves to a new CSV. (HERE I used 1 cm path length cuvette, you might have used different, so adjust accordingly)
    """
    # Molar weight (g/mol) and extinction coefficient (L·mol^-1·cm^-1)
    pigment_info = {
        "Myxoxanthophyll": {"MW": 731.01,  "E": 157700},
        "Zeaxanthin":      {"MW": 568.88,  "E": 140900},
        "Echinenone":      {"MW": 550.87,  "E": 118700},
        "B,B-carotene":    {"MW": 536.873, "E": 140400},
        "Chlorophyll-a":   {"MW": 893.51,  "E": 70020}
    }
    df = pd.read_csv(input_csv)
    required_cols = {"Pigment", "Extracted_Absorbance"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Input CSV must contain columns: {required_cols}")
    def beer_lambert_ug_mL(A, MW, E):
        # c (µg/mL) = (A * MW (g/mol) * 1000) / E
        if E == 0 or pd.isna(A):
            return None
        return (A * MW * 1000) / E
    concentrations = []
    for idx, row in df.iterrows():
        pigment = row["Pigment"]
        A = row["Extracted_Absorbance"]
        if (pigment in pigment_info) and pd.notna(A):
            MW = pigment_info[pigment]["MW"]
            E  = pigment_info[pigment]["E"]
            conc = beer_lambert_ug_mL(A, MW, E)
            concentrations.append(conc)
        else:
            concentrations.append(None)

    df["Concentration_ug_mL"] = concentrations
    df.to_csv(output_csv, index=False)
    print(f"Concentration calculation done. Results saved to '{output_csv}'.")

###############################
# MAIN: COMBINE ALL 3 STEPS  #
###############################
if __name__ == "__main__":
    print("\n========== STEP 1: Running Least-Squares Fitting Hang Tight  ==========")
    run_least_squares_fitting()
    print("\n========== STEP 2: Extract Peak Absorbances Getting There  ==========")
    extract_peak_absorbances(
        input_csv="combined_partial_contributions_optimized.csv",
        output_csv="extracted_peak_absorbances.csv"
    )
   
    print("\n========== STEP 3: Compute Concentrations Yayyy  ==========")
    compute_concentration_ug_per_mL(
        input_csv="extracted_peak_absorbances.csv",
        output_csv="with_concentrations.csv"
    )
    print("\nAll steps completed successfully, congratulations!")

