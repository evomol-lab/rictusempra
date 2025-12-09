import streamlit as st
import subprocess
import os
import shutil
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
from openbabel import pybel
from streamlit_molstar import st_molstar

# Try to import pkasolver
try:
    from pkasolver.query import calculate_microstate_pka_values
    HAS_PKASOLVER = True
    PKASOLVER_ERROR = None
except ImportError as e:
    print(f"DEBUG: pkasolver import failed: {e}")
    HAS_PKASOLVER = False
    PKASOLVER_ERROR = str(e)

# --- Helper Functions ---

def smiles_to_2d_image(smiles: str):
    """Generates a 2D image from a SMILES string using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return MolToImage(mol, size=(400, 400))
    except Exception:
        return None

def smiles_to_mol2_file(smiles: str, filename: str) -> str | None:
    """
    Generates a 3D structure, saves it as a MOL2 file, and returns the path.
    """
    try:
        mol = pybel.readstring("smi", smiles)
        mol.make3D()
        mol.write("mol2", filename, overwrite=True)
        return filename
    except Exception:
        return None

def calculate_major_microspecies_pkasolver(smiles: str, target_ph: float) -> tuple[str | None, list]:
    """
    Uses pkasolver to determine the major microspecies at a target pH.
    Returns: (best_smiles, debug_info_list)
    """
    if not HAS_PKASOLVER:
        return None, ["pkasolver library not installed."]

    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None, ["Invalid SMILES."]
        
        # Calculate pKa microstates
        # This returns a list of 'States' objects (transitions)
        states = calculate_microstate_pka_values(mol)
        
        if not states:
            # No ionizable groups found, return original
            return smiles, ["No ionizable groups identified by pkasolver."]

        # Identification Logic: "Minimal Violation Score"
        # We want the species that is consistent with the most pKa transitions at the target pH.
        
        all_smiles = set()
        for s in states:
            all_smiles.add(Chem.MolToSmiles(s.protonated_mol))
            all_smiles.add(Chem.MolToSmiles(s.deprotonated_mol))
            
        best_smiles = None
        min_violations = float('inf')
        debug_info = []

        # We also want to prefer the one closest to the pH 7 start if ties?
        # But stability is paramount.
        
        candidates = []

        for s_smiles in all_smiles:
            violations = 0
            relevant_pikas = []
            
            for state in states:
                # Directions: protonated -> deprotonated (pKa)
                # If target_ph < state.pka: Acid (Protonated) favored
                # If target_ph > state.pka: Base (Deprotonated) favored
                
                prot_smi = Chem.MolToSmiles(state.protonated_mol)
                deprot_smi = Chem.MolToSmiles(state.deprotonated_mol)
                
                is_prot = (prot_smi == s_smiles)
                is_deprot = (deprot_smi == s_smiles)
                
                if is_prot:
                    if target_ph > state.pka:
                        violations += 1 # Should be deprotonated
                        relevant_pikas.append(f"Violation: pH {target_ph} > pKa {state.pka:.2f} (should be deprotonated)")
                elif is_deprot:
                    if target_ph < state.pka:
                        violations += 1 # Should be protonated
                        relevant_pikas.append(f"Violation: pH {target_ph} < pKa {state.pka:.2f} (should be protonated)")
            
            candidates.append((violations, s_smiles))
            
            if violations < min_violations:
                min_violations = violations
                best_smiles = s_smiles
                # Store debug info for the best one found so far (or overwrite later)

        # Re-evaluate debug info for the winner
        if best_smiles:
            final_info = [f"Selected based on {min_violations} violations."]
            # Add pKa context
            for state in states:
                 prot_smi = Chem.MolToSmiles(state.protonated_mol)
                 deprot_smi = Chem.MolToSmiles(state.deprotonated_mol)
                 if prot_smi == best_smiles:
                     final_info.append(f"State is Protonated form of transition pKa {state.pka:.2f}")
                 if deprot_smi == best_smiles:
                     final_info.append(f"State is Deprotonated form of transition pKa {state.pka:.2f}")
            return best_smiles, final_info
            
        return None, ["Could not determine best species."]

    except Exception as e:
        return None, [f"Error in pkasolver: {e}"]


def run_dimorphite(smiles: str, min_ph: float, max_ph: float) -> list[str]:
    """
    Runs Dimorphite-DL and returns A LIST of all protonated SMILES strings found.
    """
    try:
        output_filename = 'protonated.smi'
        # Ensure we use the dimorphite_dl executable if in path, or python module format
        # Since we installed via pip, 'dimorphite_dl' should be in path.
        
        # If running from inside a venv where dimorphite_dl is installed:
        cmd_base = 'dimorphite_dl'
        if not shutil.which(cmd_base):
             # Fallback if not in path, try python -m
             cmd = [os.sys.executable, '-m', 'dimorphite_dl', '--ph_min', str(min_ph), '--ph_max', str(max_ph), '--output_file', output_filename, smiles]
        else:
             cmd = [cmd_base, '--ph_min', str(min_ph), '--ph_max', str(max_ph), '--output_file', output_filename, smiles]

        subprocess.run(cmd, capture_output=True, text=True, check=True)

        results = []
        if os.path.exists(output_filename):
            with open(output_filename, 'r') as f:
                for line in f:
                    if line.strip():
                        smi = line.split()[0]
                        if smi not in results:
                            results.append(smi)
            os.remove(output_filename)
            return results
    except subprocess.CalledProcessError as e:
        st.error(f"Dimorphite-DL Error: {e.stderr}")
    except Exception as e:
        # st.error(f"An error occurred: {e}")
        # Try importing directly if subprocess failed (optional fallback)
        pass
    
    if os.path.exists('protonated.smi'):
        os.remove('protonated.smi')
        
    return []

def read_file_content(filepath: str) -> str:
    """Reads and returns the content of a file."""
    with open(filepath, 'r') as f:
        return f.read()

# --- Streamlit App ---

if __name__ == "__main__":
    st.set_page_config(page_title="Molecule Protonation Tool", layout="wide")
    st.title("Rictusempra - Interactive Molecule Protonation Tool 🧪")

    # --- Sidebar for Controls ---

    with st.sidebar:
        st.image("Rictusempra.png", width="stretch")
        st.sidebar.markdown('*Rictusempra* is a web-based cheminformatics tool for interactively visualizing small molecules and calculating their most likely protonation state at a given physiological pH. \n'
                      'It provides a simple interface to generate 2D and 3D molecular structures and prepare them for further computational chemistry tasks like molecular docking or simulation.')
        st.sidebar.markdown('Please see the [documentation](https://github.com/jpmslima/rictusempra) for more information.')
        st.sidebar.markdown('Developed by the [EvoMol-Lab](github.com/evomol-lab).\n'
                        '[BioME](bioinfo.imd.ufrn.br), UFRN, Brazil.')
        st.header("Controls")
        
        smiles_input = st.text_input("Enter SMILES string:", "c1ccccc1C(=O)O")
        
        with st.form("protonation_form"):
            st.write("Set pH for Protonation")
            target_ph = st.number_input("Target pH", value=7.4, min_value=0.0, max_value=14.0, step=0.1)
            
            # Optional range for Dictorphite (hidden/optional if using pkasolver)
            use_range = st.checkbox("Define pH Range (for Dimorphite-DL)", value=False)
            if use_range:
                min_ph = st.number_input("Minimum pH", value=target_ph-0.2, min_value=0.0, max_value=14.0, step=0.1)
                max_ph = st.number_input("Maximum pH", value=target_ph+0.2, min_value=0.0, max_value=14.0, step=0.1)
            else:
                min_ph = target_ph - 0.2
                max_ph = target_ph + 0.2
                
            method = st.radio("Calculation Method", 
                              ["Advanced MicroPka (pkasolver)", "Standard (Dimorphite-DL)"],
                              index=0 if HAS_PKASOLVER else 1)
            
            if not HAS_PKASOLVER and method == "Advanced MicroPka (pkasolver)":
                st.warning(f"pkasolver is not installed. Falling back to Dimorphite-DL. Error: {PKASOLVER_ERROR}")
            
            submitted = st.form_submit_button("Calculate Protonation State")

    # --- Main Panel for Results ---

    if smiles_input:
        # --- Initial Structure Display ---
        st.header("Initial Structure", divider="rainbow")
        
        initial_mol2_path = smiles_to_mol2_file(smiles_input, "initial.mol2")
        img_initial = smiles_to_2d_image(smiles_input)

        if initial_mol2_path and img_initial:
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("2D Structure")
                st.image(img_initial, width="stretch")
            with col2:
                st.subheader("3D Structure")
                st_molstar(initial_mol2_path, key="molstar_initial")
                
                mol2_content = read_file_content(initial_mol2_path)
                st.download_button(
                    label="Download .mol2 File",
                    data=mol2_content,
                    file_name="initial_structure.mol2",
                    mime="chemical/x-mol2"
                )
        else:
            st.error("Invalid SMILES string. Please check your input.")

        # --- Protonated Structure Display (if calculation was run) ---
        if submitted:
            st.header("Protonated Structure(s)", divider="rainbow")
            
            results = [] # List of (smiles, label)
            
            if method == "Advanced MicroPka (pkasolver)" and HAS_PKASOLVER:
                with st.spinner("Running pkasolver (Micro-pKa)..."):
                    best_smi, info = calculate_major_microspecies_pkasolver(smiles_input, target_ph)
                    if best_smi:
                        results.append((best_smi, "Major Microspecies (pkasolver)"))
                        with st.expander("See pkasolver details"):
                            for i in info:
                                st.write(f"- {i}")
                    else:
                         st.error("pkasolver failed to return a structure.")
                         for i in info:
                                st.write(f"- {i}")

            else: # Dimorphite
                with st.spinner("Running Dimorphite-DL..."):
                    smi_list = run_dimorphite(smiles_input, min_ph, max_ph)
                    if smi_list:
                        for i, smi in enumerate(smi_list):
                            results.append((smi, f"Variant {i+1} (Dimorphite-DL)"))
                    else:
                        st.warning("Dimorphite-DL found no protonation states in this range.")

            # --- Display Results ---
            if results:
                st.success(f"Found {len(results)} relevant structure(s).")
                
                # If multiple, use tabs
                if len(results) > 1:
                    tabs = st.tabs([label for _, label in results])
                    for i, tab in enumerate(tabs):
                        smi, label = results[i]
                        with tab:
                            st.code(smi, language="smiles")
                            mol2_path = smiles_to_mol2_file(smi, f"protonated_{i}.mol2")
                            img = smiles_to_2d_image(smi)
                            
                            if mol2_path and img:
                                c3, c4 = st.columns(2)
                                with c3:
                                    st.image(img, width="stretch")
                                with c4:
                                    st_molstar(mol2_path, key=f"molstar_res_{i}")
                                    content = read_file_content(mol2_path)
                                    st.download_button("Download", content, f"structure_{i}.mol2", "chemical/x-mol2", key=f"dl_{i}")
                else:
                    # Single result
                    smi, label = results[0]
                    st.subheader(label)
                    st.code(smi, language="smiles")
                    mol2_path = smiles_to_mol2_file(smi, f"protonated.mol2")
                    img = smiles_to_2d_image(smi)
                    
                    if mol2_path and img:
                        c3, c4 = st.columns(2)
                        with c3:
                            st.image(img, width="stretch")
                        with c4:
                            st_molstar(mol2_path, key=f"molstar_res_single")
                            content = read_file_content(mol2_path)
                            st.download_button("Download", content, "protonated_structure.mol2", "chemical/x-mol2")
