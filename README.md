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

---

## Core Features

- **SMILES Input**: Accepts a SMILES string to define the initial molecule.

- **2D & 3D Visualization**: Instantly renders 2D chemical diagrams (via RDKit) and interactive 3D structures (via Open Babel & streamlit-molstar).

- **Protonation State Calculation**:
    - **Advanced**: Uses `pkasolver` to identify the major microspecies based on specific pKa transitions.
    - **Standard**: Uses `Dimorphite-DL` to enumerate highly probable protonation states within a user-defined pH range.

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

##  <a name='Disclaimer'></a>Disclaimer

The developer team used generative AI tools for the following tasks:

- Code revision and optimization.

- Elaborate documentation topic structure.

- Review english language.

----

<div align="center">
  <img src="EvoMol_v2.png" alt="Rictusempra Logo" width="400">
</div>