<div align="center">
  <img src="Rictusempra.png" alt="Rictusempra Logo" width="300">
</div>

# ***Rictusempra***: Interactive Molecule Protonation Tool 🧪

Rictusempra is a web-based cheminformatics tool for interactively visualizing small molecules and calculating their most likely protonation state at a given physiological pH. It provides a simple interface to generate 2D and 3D molecular structures and prepare them for further computational chemistry tasks like molecular docking or simulation.

<div align="center">
  <img src="EvoMol-logo.png" alt="Rictusempra Logo" width="100">
</div>

Developed by the [EvoMol-Lab](github.com/evomol-lab), [BioME](bioinfo.imd.ufrn.br), UFRN, Brazil.

> **Note:** The application's user interface (labels, charts, tables) is entirely in **English**. This documentation is also in English; only commit messages and pull request discussions for this project may be written in Portuguese.

---

## Highlights & Improvements 🚀

- **Advanced Micro-pKa Prediction**: Now integrates **pkasolver**, a Graph Neural Network (GNN) model, to scientifically predict micro-pKa values and determining the major microspecies with high accuracy.
- **Enhanced Fallback**: Improved integration with **Dimorphite-DL** to display *all* plausible protonation states (microspecies) when using the standard method.
- **Interactive Methodology**: Choose between "Advanced (pkasolver)" and "Standard (Dimorphite-DL)" calculation modes directly from the UI.
- **Tautomer/Microspecies Ratio vs. pH**: For each ionizable group identified by pkasolver, plots the protonated/deprotonated ratio across a full pH range (Henderson-Hasselbalch), with the group's pKa and the chosen target pH marked on the curve, plus a summary table with the exact percentages at the target pH.
- **Microspecies Distribution at the Target pH**: Splits the molecule into its sequential microspecies (from the most protonated to the most deprotonated form, as defined by the pKa ladder) and shows the exact percentage of each one at the chosen pH — the percentages always add up to 100%.
- **Microspecies Distribution Across the Full pH Range**: A classic species-distribution diagram — a stacked area chart sweeping pH 0–14 in 0.5 steps, one color per microspecies, always summing to 100% at every pH — plus a downloadable CSV of the underlying data.
- **Microspecies Gallery with Structure Downloads**: Below the microspecies table, small 2D structure images of every microspecies on the pKa ladder are shown side by side, numbered to match the table and organized in rows. The major microspecies at the target pH is highlighted with a gold border. Each structure's 3D `.sdf` can be downloaded individually, and all of them can also be downloaded together as a single multi-structure `.sdf` file.

---

## Core Features

- **SMILES Input**: Accepts a SMILES string to define the initial molecule.

- **2D & 3D Visualization**: Instantly renders 2D chemical diagrams and interactive 3D structures (via RDKit & streamlit-molstar).

- **Protonation State Calculation**:
    - **Advanced**: Uses `pkasolver` to identify the major microspecies based on specific pKa transitions.
    - **Standard**: Uses `Dimorphite-DL` to enumerate highly probable protonation states within a user-defined pH range.

- **Tautomer/Microspecies Ratio vs. pH and pKa**: When using the Advanced (pkasolver) method, shows an interactive chart with the protonated/deprotonated ratio of each ionizable group as a function of pH (Henderson-Hasselbalch), plus a table with the exact percentages at the selected target pH.

- **Microspecies Distribution at the Target pH**: Also under the Advanced (pkasolver) method, a second chart shows how the molecule's total population is split across its sequential microspecies (most protonated → most deprotonated) at the chosen pH, chaining Henderson-Hasselbalch along the identified pKa ladder — the percentages shown always sum to 100%, unlike the per-site ratio above, which considers each ionizable group in isolation.

- **Microspecies Distribution Across the Full pH Range**: A third chart sweeps the same microspecies breakdown across pH 0–14 in steps of 0.5, rendered as a stacked area chart with one color per microspecies (stack always sums to 100% at every pH), with the underlying pH-sweep table available as a CSV download.

- **Microspecies Gallery with Structure Downloads**: Under the microspecies table, every microspecies on the pKa ladder — from the most protonated to the most deprotonated form — gets its own small 2D structure image, numbered to match the table and laid out in rows. The major microspecies at the target pH is visually highlighted (gold border, and marked with ⭐ in both the table and the gallery) so it stays easy to spot among the rest. Each structure's 3D `.sdf` can be downloaded individually, or all of them at once as a single multi-structure `.sdf` file.

- **Side-by-Side Comparison**: Displays the initial and protonated structures next to each other for easy comparison.

- **Structure Download**: Allows users to download the generated 3D structures in `.sdf` format, for the initial structure, each protonation-state result, and every individual microspecies (or all microspecies at once).

---

## Technology Stack

- **Frontend**: Streamlit

- **2D Structure Rendering**: RDKit

- **3D Structure Generation**: RDKit (ETKDG embedding + MMFF/UFF force field optimization)

- **3D Structure Visualization**: streamlit-molstar

- **Protonation Calculation**:
    - **pkasolver** (Graph Neural Networks)
    - **Dimorphite-DL** (Rule-based)

- **Ratio-vs-pH Charts**: Pandas & Altair

---

## Installation and Setup

It is **highly recommended** to use Conda for installation, as it handles the complex dependencies of RDKit smoothly.

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
conda install -c conda-forge rdkit dimorphite-dl
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
    - You can download the `.sdf` files for any generated structure.

6. **Check the Tautomer/Microspecies Ratio (Advanced method only)**: Below the pkasolver results, the "Ratio of each tautomer as a function of pH and pKa" section shows a chart with one curve per ionizable group — its protonated/deprotonated ratio across pH 0–14, with its pKa and your target pH marked — plus a table with the exact percentages at the target pH.

7. **Check the Microspecies Distribution (Advanced method only)**: Right below that, the "Amount of each microspecies at pH X" section shows a bar chart with the percentage of each sequential microspecies (most protonated → most deprotonated) at your target pH, plus a table and the sum of the percentages (always 100%). The row for the major microspecies at the target pH is highlighted in gold.

8. **View and Download Every Microspecies Structure (Advanced method only)**: Below the table, the "Microspecies Structures" section shows a small 2D structure image for each microspecies, numbered to match the table and organized in rows. The major microspecies is outlined in gold and marked with ⭐. Use the "Download .sdf" button under any structure to get its individual 3D file, or the "Download all microspecies (.sdf, multi-structure)" button above the gallery to get every microspecies in one combined `.sdf` file.

9. **Explore the Full pH Sweep (Advanced method only)**: The "Microspecies distribution across the pH range" section shows the same breakdown as a stacked area chart swept over pH 0–14 (steps of 0.5), one color per microspecies. Use the "Download pH sweep data (.csv)" button to get the full table (pH, microspecies, description, fraction, percentage) for further analysis.

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

## License

Rictusempra's own source code is released under the **MIT License** — see [`LICENSE`](LICENSE).

### Third-party licenses

Rictusempra is built on top of several open-source packages, each under its own license:

| Package | License |
|---|---|
| [RDKit](https://www.rdkit.org) | BSD-3-Clause |
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

All of Rictusempra's dependencies are permissive (MIT/BSD/Apache-2.0) or weak-copyleft (LGPL-3.0, used by `cairosvg`, which permits dynamic use without requiring the rest of the application to be LGPL). Earlier versions of this tool used **Open Babel (GPL-2.0)** for 3D structure generation; that dependency has been replaced with RDKit's own ETKDG embedding + MMFF/UFF optimization, so no strong-copyleft (GPL) component remains anywhere in the stack.

##  <a name='Disclaimer'></a>Disclaimer

The developer team used generative AI tools for the following tasks:

- Code revision and optimization.

- Elaborate documentation topic structure.

- Review english language.

----

<div align="center">
  <img src="EvoMol_v2.png" alt="Rictusempra Logo" width="400">
</div>