<div align="center">
  <img src="Rictusempra.png" alt="Rictusempra Logo" width="300">
</div>

# ***Rictusempra***: Interactive Molecule Protonation Tool 🧪

Rictusempra is a web-based cheminformatics tool for interactively visualizing small molecules and calculating their most likely protonation state at a given physiological pH. It provides a simple interface to generate 2D and 3D molecular structures and prepare them for further computational chemistry tasks like molecular docking or simulation.

<div align="center">
  <img src="EvoMol-logo.png" alt="Rictusempra Logo" width="100">
</div>

Developed by the [EvoMol-Lab](github.com/evomol-lab), [BioME](bioinfo.imd.ufrn.br), UFRN, Brazil.

---

## Highlights & Improvements 🚀

- **Advanced Micro-pKa Prediction**: Now integrates **pkasolver**, a Graph Neural Network (GNN) model, to scientifically predict micro-pKa values and determining the major microspecies with high accuracy.
- **Enhanced Fallback**: Improved integration with **Dimorphite-DL** to display *all* plausible protonation states (microspecies) when using the standard method.
- **Interactive Methodology**: Choose between "Advanced (pkasolver)" and "Standard (Dimorphite-DL)" calculation modes directly from the UI.
- **Tautomer/Microspecies Ratio vs. pH**: For each ionizable group identified by pkasolver, plots the protonated/deprotonated ratio across a full pH range (Henderson-Hasselbalch), with the group's pKa and the chosen target pH marked on the curve, plus a summary table with the exact percentages at the target pH.

---

## Core Features

- **SMILES Input**: Accepts a SMILES string to define the initial molecule.

- **2D & 3D Visualization**: Instantly renders 2D chemical diagrams (via RDKit) and interactive 3D structures (via Open Babel & streamlit-molstar).

- **Protonation State Calculation**:
    - **Advanced**: Uses `pkasolver` to identify the major microspecies based on specific pKa transitions.
    - **Standard**: Uses `Dimorphite-DL` to enumerate highly probable protonation states within a user-defined pH range.

- **Tautomer/Microspecies Ratio vs. pH and pKa**: When using the Advanced (pkasolver) method, shows an interactive chart with the protonated/deprotonated ratio of each ionizable group as a function of pH (Henderson-Hasselbalch), plus a table with the exact percentages at the selected target pH.

- **Side-by-Side Comparison**: Displays the initial and protonated structures next to each other for easy comparison.

- **Structure Download**: Allows users to download the generated 3D structures in `.mol2` format.

---

## Technology Stack

- **Frontend**: Streamlit

- **2D Structure Rendering**: RDKit

- **3D Structure Generation**: Open Babel

- **3D Structure Visualization**: streamlit-molstar

- **Protonation Calculation**:
    - **pkasolver** (Graph Neural Networks)
    - **Dimorphite-DL** (Rule-based)

- **Ratio-vs-pH Charts**: Pandas & Altair

---

## Installation and Setup

It is **highly recommended** to use Conda for installation, as it handles the complex dependencies of RDKit and Open Babel smoothly.

### Step 1: Clone the Repository

```
git clone https://github.com/evomol-lab/rictusempra.git
cd rictusempra
```

### Step 2: Create Conda Environment

It's best practice to create a dedicated environment for the tool.

```
# Create and activate the conda environment
conda create -n rictusempra python=3.9
conda activate rictusempra

# Install packages from conda-forge
conda install -c conda-forge rdkit openbabel dimorphite-dl
```

### Step 3: Install Pip Dependencies

Install the remaining Python packages using pip.

```
pip install streamlit streamlit-molstar torch torch-geometric
```

*Note: For `pkasolver`, you may need to install it from source or check compatibility with your specific environment.*

---

## Usage

Once the environment is set up, you can run the Streamlit application from your terminal.

```
streamlit run rictusempra.py
```

A new tab will open in your web browser with the application running.

### How to Use the Tool:

1. **Enter SMILES**: In the sidebar on the left, enter the SMILES string of the molecule you want to analyze. The initial 2D and 3D structures will appear on the main panel.

2. **Set pH**: Enter the target pH (default 7.4).
    - If using **Dimorphite-DL**, you can optionally define a min/max range.

3. **Select Method**: Choose between:
    - **Advanced MicroPka (pkasolver)**: Best for finding the single most stable major microspecies.
    - **Standard (Dimorphite-DL)**: Best for finding a range of possible states.

4. **Calculate**: Click the "Calculate Protonation State" button.

5. **View & Download**: The results for the protonated molecule(s) will appear below.
    - If multiple valid states are found, they will be shown in tabs.
    - You can download the `.mol2` files for any generated structure.

6. **Check the Tautomer/Microspecies Ratio (Advanced method only)**: Below the pkasolver results, the "Razão de cada tautômero em função do pH e do pKa" section shows a chart with one curve per ionizable group — its protonated/deprotonated ratio across pH 0–14, with its pKa and your target pH marked — plus a table with the exact percentages at the target pH.

---

## Project File Structure

For the application to work correctly, your project folder should be organized as follows:

```
/your-project-folder
|-- rictusempra.py       # The main Streamlit app script
|-- rictusempra.png      # Your logo file
|-- requirements.txt   # Python dependencies
|-- packages.txt       # System dependencies
|-- README.md          # This documentation file
|-- pkasolver/         # (Optional) Local installation of pkasolver if used
```

---

## References

If you use this tool in your research, please cite the underlying open-source packages that make it possible:

- **pkasolver**:
  Mayr, F., Wieder, O., Wieder, M., & Langer, T. (2022). Improving Small Molecule pKa Prediction Using Transfer Learning with Graph Neural Networks. *bioRxiv*. [https://doi.org/10.1101/2022.01.20.476787](https://www.biorxiv.org/content/10.1101/2022.01.20.476787v1)

- **Dimorphite-DL**:
  Ropp, P. J., Kaminsky, J. C., Yablonski, S., & Durrant, J. D. (2019). Dimorphite-DL: An open-source program for enumerating the ionization states of drug-like small molecules. *Journal of Cheminformatics*, *11*(1), 51. [https://doi.org/10.1186/s13321-019-0371-5](https://doi.org/10.1186/s13321-019-0371-5).

- **RDKit**:
  RDKit: Open-Source Cheminformatics Software. (n.d.). Retrieved August 24, 2025, from [http://www.rdkit.org](http://www.rdkit.org)

- **Open Babel**:
  O'Boyle, N. M., Banck, M., James, C. A., Morley, C., Vandermeersch, T., & Hutchison, G. R. (2011). Open Babel: An open chemical toolbox. *Journal of Cheminformatics*, *3*(1), 33. [https://doi.org/10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33)

## License

Rictusempra's own source code is released under the **MIT License** — see [`LICENSE`](LICENSE).

### Third-party licenses

Rictusempra is built on top of several open-source packages, each under its own license:

| Package | License |
|---|---|
| [RDKit](https://www.rdkit.org) | BSD-3-Clause |
| [Open Babel](https://openbabel.org) (`openbabel-wheel`) | **GPL-2.0** |
| [Dimorphite-DL](https://github.com/durrantlab/dimorphite_dl) | Apache-2.0 |
| [pkasolver](https://github.com/mayrf/pkasolver) (vendored in this repo) | MIT |
| [Streamlit](https://streamlit.io) | Apache-2.0 |
| [streamlit-molstar](https://github.com/pragmatic-streamlit/streamlit-molstar) | OSI-approved (see upstream repo) |
| [PyTorch](https://pytorch.org) | BSD-3-Clause |
| [PyTorch Geometric](https://pytorch-geometric.readthedocs.io) | MIT |
| [CairoSVG](https://cairosvg.org) | LGPL-3.0-or-later |
| [svgutils](https://svgutils.readthedocs.io) | MIT |
| [pandas](https://pandas.pydata.org) | BSD-3-Clause |
| [Vega-Altair](https://altair-viz.github.io) | BSD-3-Clause |
| [NumPy](https://numpy.org) | BSD-3-Clause |

**Note on Open Babel (GPL-2.0):** Rictusempra imports Open Babel's Python bindings (`openbabel.pybel`) directly in its source code. Because of this, if you redistribute Rictusempra as a combined, runnable application (rather than just reusing its MIT-licensed original code on its own), that combined distribution is subject to the terms of the GPL-2.0 for the Open Babel component. If this matters for your use case (e.g., embedding Rictusempra in proprietary software), consult the GPL-2.0 terms or consider swapping Open Babel for a non-copyleft alternative.

##  <a name='Disclaimer'></a>Disclaimer

The developer team used generative AI tools for the following tasks:

- Code revision and optimization.

- Elaborate documentation topic structure.

- Review english language.

----

<div align="center">
  <img src="EvoMol_v2.png" alt="Rictusempra Logo" width="400">
</div>