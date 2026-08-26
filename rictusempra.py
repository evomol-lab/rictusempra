import streamlit as st
import subprocess
import os
import shutil
import io
import numpy as np
import pandas as pd
import altair as alt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolToImage
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

def _embed_and_optimize_3d(mol: Chem.Mol) -> Chem.Mol | None:
    """
    Embeds 3D coordinates (ETKDG) and optimizes the geometry (MMFF94,
    falling back to UFF) for an RDKit Mol that already has explicit
    hydrogens. Mutates and returns `mol`, or returns None if embedding
    fails even after the random-coordinates fallback.

    Uses RDKit instead of Open Babel so the whole toolchain stays under
    permissive (non-copyleft) licenses.
    """
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    if AllChem.EmbedMolecule(mol, params) != 0:
        # Retry with random coordinates as a fallback for tricky molecules
        if AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=0xF00D) != 0:
            return None

    # Optimize geometry; MMFF94 first, falling back to UFF if unsupported
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception:
            pass

    return mol

def smiles_to_3d_file(smiles: str, filename: str) -> str | None:
    """
    Generates a 3D structure using RDKit (ETKDG embedding + MMFF/UFF force
    field optimization) and saves it as an SDF file. Returns the path.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        if _embed_and_optimize_3d(mol) is None:
            return None

        writer = Chem.SDWriter(filename)
        writer.write(mol)
        writer.close()
        return filename
    except Exception:
        return None

def get_smiles_from_pubchem(cid: int) -> str | None:
    """Fetches the SMILES string for a given PubChem CID."""
    import urllib.request
    import json
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
    try:
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())
            props = data['PropertyTable']['Properties'][0]
            return props.get('CanonicalSMILES') or props.get('IsomericSMILES') or props.get('SMILES') or props.get('ConnectivitySMILES')
    except Exception:
        return None

def calculate_major_microspecies_pkasolver(smiles: str, target_ph: float) -> tuple[str | None, list, list]:
    """
    Uses pkasolver to determine the major microspecies at a target pH.
    Returns: (best_smiles, debug_info_list, states)
    `states` is the list of pkasolver 'States' objects (one per ionizable
    group/tautomeric transition identified), used to plot the ratio of each
    species as a function of pH.
    """
    if not HAS_PKASOLVER:
        return None, ["pkasolver library not installed."], []

    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None, ["Invalid SMILES."], []

        # Calculate pKa microstates
        # This returns a list of 'States' objects (transitions)
        states = calculate_microstate_pka_values(mol)

        if not states:
            # No ionizable groups found, return original
            return smiles, ["No ionizable groups identified by pkasolver."], []

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
            return best_smiles, final_info, states

        return None, ["Could not determine best species."], states

    except Exception as e:
        return None, [f"Error in pkasolver: {e}"], []


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

def protonated_fraction(ph, pka):
    """
    Fraction of the protonated form of an ionizable group at a given pH,
    using the Henderson-Hasselbalch equation:

        fraction_protonated = 1 / (1 + 10 ** (pH - pKa))

    `ph` can be a scalar or a numpy array.
    """
    return 1.0 / (1.0 + 10.0 ** (np.asarray(ph, dtype=float) - pka))

def build_ph_distribution_df(states: list, ph_min: float = 0.0, ph_max: float = 14.0, num_points: int = 281) -> pd.DataFrame:
    """
    Builds a long-format DataFrame with the ratio (fraction of the population)
    of the protonated and deprotonated form of each ionizable group/tautomer
    identified by pkasolver, across a range of pH values.

    Columns: pH, Site (site label including its pKa), pKa, Form
    (Protonated/Deprotonated), Fraction (0-1).
    """
    ph_values = np.linspace(ph_min, ph_max, num_points)
    rows = []
    for idx, state in enumerate(states, start=1):
        pka = float(state.pka)
        frac_prot = protonated_fraction(ph_values, pka)
        frac_deprot = 1.0 - frac_prot
        site_label = f"Tautomer/site {idx} (pKa {pka:.2f})"
        for ph, fp, fd in zip(ph_values, frac_prot, frac_deprot):
            rows.append({"pH": float(ph), "Site": site_label, "pKa": pka, "Form": "Protonated", "Fraction": float(fp)})
            rows.append({"pH": float(ph), "Site": site_label, "pKa": pka, "Form": "Deprotonated", "Fraction": float(fd)})
    return pd.DataFrame(rows)

def summarize_fraction_at_ph(states: list, ph: float) -> pd.DataFrame:
    """
    Summarizes, for the chosen target pH, the ratio (%) of the protonated and
    deprotonated form of each ionizable group/tautomer identified by pkasolver.
    """
    rows = []
    for idx, state in enumerate(states, start=1):
        pka = float(state.pka)
        frac_prot = float(protonated_fraction(ph, pka))
        frac_deprot = 1.0 - frac_prot
        rows.append({
            "Tautomer/Site": idx,
            "pKa": round(pka, 2),
            "% Protonated": round(frac_prot * 100, 1),
            "% Deprotonated": round(frac_deprot * 100, 1),
        })
    return pd.DataFrame(rows)

def build_macrospecies_distribution_at_ph(states: list, ph: float) -> pd.DataFrame:
    """
    Computes the population distribution across the molecule's sequential
    macrospecies (overall protonation states) at a given pH.

    The `n` pKa transitions identified by pkasolver define a ladder of
    `n + 1` macrospecies, from the most protonated (index 0, before any of
    the identified groups has lost a proton) to the most deprotonated
    (index n, after all of them have). Each step between consecutive
    macrospecies follows Henderson-Hasselbalch, so the fraction of
    macrospecies `i` relative to macrospecies `i-1` is `10 ** (pH - pKa_i)`.
    Chaining these ratios and normalizing gives a distribution that always
    sums to 1 (100%) — unlike the per-site ratios in
    `build_ph_distribution_df`, which treat each site independently.
    """
    pkas = sorted(float(s.pka) for s in states)
    n = len(pkas)

    # log10 of the (unnormalized) weight of each macrospecies, built
    # cumulatively to stay numerically stable regardless of how many
    # ionizable groups are chained together.
    log_weights = [0.0]
    for pka in pkas:
        log_weights.append(log_weights[-1] + (ph - pka))

    max_log = max(log_weights)
    weights = [10.0 ** (lw - max_log) for lw in log_weights]
    total = sum(weights)
    fractions = [w / total for w in weights]

    rows = []
    for i, frac in enumerate(fractions):
        if i == 0:
            desc = "most protonated"
        elif i == n:
            desc = "most deprotonated"
        else:
            desc = f"{i} proton removed" if i == 1 else f"{i} protons removed"
        rows.append({
            "Microspecies": f"Microspecies {i + 1}",
            "Description": desc,
            "Protons removed": i,
            "Fraction": frac,
            "Percentage": frac * 100.0,
        })
    return pd.DataFrame(rows)

def plot_macrospecies_distribution(df: pd.DataFrame, target_ph: float) -> alt.LayerChart:
    """
    Bar chart of the population distribution across macrospecies at the
    target pH (bars sum to 100%), with the percentage labeled on each bar.
    """
    order = df["Microspecies"].tolist()
    bars = alt.Chart(df).mark_bar().encode(
        x=alt.X("Microspecies:N", sort=order, title="Microspecies (most protonated → most deprotonated)"),
        y=alt.Y("Percentage:Q", title="% of total population", scale=alt.Scale(domain=[0, 100])),
        tooltip=["Microspecies", "Description", alt.Tooltip("Percentage:Q", format=".2f")],
    )
    labels = bars.mark_text(dy=-8).encode(text=alt.Text("Percentage:Q", format=".1f"))
    return (bars + labels).properties(height=320, title=f"Microspecies distribution at pH {target_ph:.2f}")

def build_macrospecies_distribution_sweep(states: list, ph_min: float = 0.0, ph_max: float = 14.0, step: float = 0.5) -> pd.DataFrame:
    """
    Sweeps pH from `ph_min` to `ph_max` (inclusive) in increments of `step`,
    computing the macrospecies population distribution (see
    `build_macrospecies_distribution_at_ph`) at each point. At every pH value
    the percentages across microspecies still add up to 100%.
    """
    num_steps = int(round((ph_max - ph_min) / step))
    rows = []
    for i in range(num_steps + 1):
        ph = round(ph_min + i * step, 10)
        df = build_macrospecies_distribution_at_ph(states, ph)
        df.insert(0, "pH", round(ph, 2))
        rows.append(df)
    return pd.concat(rows, ignore_index=True)

def plot_macrospecies_sweep(df: pd.DataFrame, target_ph: float) -> alt.LayerChart:
    """
    Stacked area chart of the macrospecies population distribution across a
    full pH sweep, one color per microspecies. At every pH the stacked
    percentages add up to 100%. A dashed vertical line marks the currently
    selected target pH.
    """
    order = sorted(df["Microspecies"].unique(), key=lambda m: int(m.rsplit(" ", 1)[-1]))
    area = alt.Chart(df).mark_area().encode(
        x=alt.X("pH:Q", title="pH"),
        y=alt.Y("Percentage:Q", title="% of total population", stack="zero", scale=alt.Scale(domain=[0, 100])),
        color=alt.Color("Microspecies:N", sort=order, title="Microspecies"),
        order=alt.Order("Protons removed:Q"),
        tooltip=["Microspecies", "Description", alt.Tooltip("pH:Q", format=".2f"), alt.Tooltip("Percentage:Q", format=".2f")],
    )
    ph_rule = alt.Chart(pd.DataFrame({"pH": [target_ph]})).mark_rule(color="red", strokeDash=[4, 4]).encode(x="pH:Q")
    return (area + ph_rule).properties(height=350, title="Microspecies distribution across the pH range")

def plot_ph_distribution(df: pd.DataFrame, target_ph: float) -> alt.LayerChart:
    """
    Builds an Altair chart with one curve per ionizable group/tautomer showing
    how its protonated/deprotonated ratio changes with pH, plus reference
    lines marking each pKa and the currently selected target pH.
    """
    base = alt.Chart(df).mark_line().encode(
        x=alt.X("pH:Q", title="pH"),
        y=alt.Y("Fraction:Q", title="Ratio (fraction of the population)", scale=alt.Scale(domain=[0, 1])),
        color=alt.Color("Form:N", title="Form"),
        strokeDash=alt.StrokeDash("Site:N", title="Tautomer / ionizable site"),
        tooltip=["Site", "Form", alt.Tooltip("pH:Q", format=".2f"), alt.Tooltip("Fraction:Q", format=".2%"), alt.Tooltip("pKa:Q", format=".2f")],
    )

    pka_df = df[["Site", "pKa"]].drop_duplicates()
    pka_rule = alt.Chart(pka_df).mark_rule(strokeDash=[2, 2], color="gray").encode(
        x="pKa:Q",
        tooltip=["Site", alt.Tooltip("pKa:Q", format=".2f")],
    )

    ph_rule = alt.Chart(pd.DataFrame({"pH": [target_ph]})).mark_rule(
        color="red", strokeDash=[4, 4]
    ).encode(x="pH:Q")

    return (base + pka_rule + ph_rule).properties(height=350)

def build_microspecies_ladder(states: list) -> list[dict]:
    """
    Reconstructs the structure of every microspecies on the sequential pKa
    ladder (see `build_macrospecies_distribution_at_ph`), from the most
    protonated to the most deprotonated form.

    Sorted by ascending pKa, `state.protonated_mol` is the species favored
    below that transition's pKa and `state.deprotonated_mol` the species
    favored above it, with each transition's deprotonated form equal to the
    next transition's protonated form. Chaining them therefore reconstructs
    all `n + 1` rungs of the ladder from the `n` transitions.

    Returns a list of dicts (same order/indexing as
    `build_macrospecies_distribution_at_ph`): {"index" (1-based, matching
    the "Microspecies N" labels), "smiles", "protons_removed"}.
    """
    sorted_states = sorted(states, key=lambda s: s.pka)
    rungs = [Chem.MolToSmiles(sorted_states[0].protonated_mol)]
    for state in sorted_states:
        rungs.append(Chem.MolToSmiles(state.deprotonated_mol))

    return [
        {"index": i + 1, "smiles": smi, "protons_removed": i}
        for i, smi in enumerate(rungs)
    ]

def find_major_microspecies_index(ladder: list[dict], major_smiles: str | None) -> int | None:
    """Returns the 1-based ladder index whose SMILES matches `major_smiles`, or None."""
    if not major_smiles:
        return None
    for item in ladder:
        if item["smiles"] == major_smiles:
            return item["index"]
    return None

def build_microspecies_sdf_content(ladder: list[dict]) -> str | None:
    """
    Generates a 3D structure (same RDKit ETKDG + MMFF/UFF pipeline as
    `smiles_to_3d_file`) for every microspecies in `ladder` and writes them
    into a single multi-record .sdf string, tagging each record with its
    microspecies index, protons removed, and SMILES as SDF properties.
    """
    buffer = io.StringIO()
    writer = Chem.SDWriter(buffer)
    wrote_any = False
    for item in ladder:
        mol = Chem.MolFromSmiles(item["smiles"])
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        if _embed_and_optimize_3d(mol) is None:
            continue
        mol.SetProp("_Name", f"Microspecies_{item['index']}")
        mol.SetProp("Microspecies", str(item["index"]))
        mol.SetProp("ProtonsRemoved", str(item["protons_removed"]))
        mol.SetProp("SMILES", item["smiles"])
        writer.write(mol)
        wrote_any = True
    writer.close()
    return buffer.getvalue() if wrote_any else None

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
        
        input_type = st.radio("Input Method:", ["SMILES", "PubChem CID"])
        
        if input_type == "SMILES":
            smiles_input = st.text_input("Enter SMILES string:", "C1=C(NC=N1)CC(C(=O)O)N")
        else:
            cid_input = st.number_input("Enter PubChem CID:", value=773, step=1, min_value=1)
            smiles_input = get_smiles_from_pubchem(cid_input)
            if smiles_input:
                st.caption(f"Loaded SMILES: {smiles_input}")
            else:
                st.error("Failed to load SMILES from PubChem")
        
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
        
        initial_structure_path = smiles_to_3d_file(smiles_input, "initial.sdf")
        img_initial = smiles_to_2d_image(smiles_input)

        if initial_structure_path and img_initial:
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("2D Structure")
                st.image(img_initial, width="stretch")
            with col2:
                st.subheader("3D Structure")
                st_molstar(initial_structure_path, key="molstar_initial")

                sdf_content = read_file_content(initial_structure_path)
                st.download_button(
                    label="Download .sdf File",
                    data=sdf_content,
                    file_name="initial_structure.sdf",
                    mime="chemical/x-mdl-sdfile"
                )
        else:
            st.error("Invalid SMILES string. Please check your input.")

        # --- Protonated Structure Display (if calculation was run) ---
        if submitted:
            st.header("Protonated Structure(s)", divider="rainbow")
            
            results = [] # List of (smiles, label)
            
            if method == "Advanced MicroPka (pkasolver)" and HAS_PKASOLVER:
                with st.spinner("Running pkasolver (Micro-pKa)..."):
                    best_smi, info, states = calculate_major_microspecies_pkasolver(smiles_input, target_ph)
                    if best_smi:
                        results.append((best_smi, "Major Microspecies (pkasolver)"))
                        with st.expander("See pkasolver details"):
                            for i in info:
                                st.write(f"- {i}")
                    else:
                         st.error("pkasolver failed to return a structure.")
                         for i in info:
                                st.write(f"- {i}")

                    # --- Ratio of each tautomer/microspecies vs pH and pKa ---
                    if states:
                        st.subheader("Ratio of each tautomer as a function of pH and pKa")
                        st.caption(
                            "For each ionizable group (tautomer/microspecies) identified by "
                            "pkasolver, the curve shows the ratio (fraction of the population) of "
                            "the protonated and deprotonated forms across pH, calculated with the "
                            "Henderson-Hasselbalch equation. The gray dashed line marks each "
                            "group's pKa; the red dashed line marks the selected target pH."
                        )
                        dist_df = build_ph_distribution_df(states, ph_min=0.0, ph_max=14.0)
                        st.altair_chart(plot_ph_distribution(dist_df, target_ph), use_container_width=True)
                        st.dataframe(
                            summarize_fraction_at_ph(states, target_ph),
                            width="stretch",
                            hide_index=True,
                        )

                        # --- Population distribution across macrospecies at the target pH (sums to 100%) ---
                        st.subheader(f"Amount of each microspecies at pH {target_ph:.2f}")
                        st.caption(
                            "Unlike the chart above (which treats each site independently), here "
                            "the molecule is split into its sequential microspecies defined by the "
                            "identified pKa values — from the most protonated to the most "
                            "deprotonated form — and the fraction of each one at the chosen pH is "
                            "calculated by chaining Henderson-Hasselbalch along the pKa ladder. "
                            "The fractions always add up to 1 (100%)."
                        )
                        macro_df = build_macrospecies_distribution_at_ph(states, target_ph)
                        st.altair_chart(plot_macrospecies_distribution(macro_df, target_ph), use_container_width=True)

                        ladder = build_microspecies_ladder(states)
                        major_index = find_major_microspecies_index(ladder, best_smi)

                        table_df = macro_df[["Microspecies", "Description", "Percentage"]].round({"Percentage": 2})
                        if major_index is not None:
                            def _highlight_major(row):
                                is_major = int(row["Microspecies"].rsplit(" ", 1)[-1]) == major_index
                                style = "background-color: rgba(255, 215, 0, 0.25); font-weight: bold"
                                return [style if is_major else "" for _ in row]
                            st.dataframe(table_df.style.apply(_highlight_major, axis=1), width="stretch", hide_index=True)
                        else:
                            st.dataframe(table_df, width="stretch", hide_index=True)
                        st.caption(f"Sum of fractions: {macro_df['Percentage'].sum():.2f}%")

                        # --- Individual microspecies structures: view & download .sdf ---
                        st.subheader("Microspecies Structures")
                        st.caption(
                            "2D structure of every microspecies on the pKa ladder above, numbered to "
                            "match the table (Microspecies 1 = most protonated → last = most "
                            "deprotonated). The major microspecies at the target pH is highlighted in "
                            "gold. Each structure can be downloaded as a 3D .sdf file individually, or "
                            "all together below."
                        )
                        if major_index is not None:
                            st.markdown(f"⭐ **Major microspecies at pH {target_ph:.2f}: Microspecies {major_index}**")

                        combined_sdf = build_microspecies_sdf_content(ladder)
                        if combined_sdf:
                            st.download_button(
                                label="Download all microspecies (.sdf, multi-structure)",
                                data=combined_sdf,
                                file_name="all_microspecies.sdf",
                                mime="chemical/x-mdl-sdfile",
                            )

                        cols_per_row = 4
                        for row_start in range(0, len(ladder), cols_per_row):
                            row_items = ladder[row_start:row_start + cols_per_row]
                            row_cols = st.columns(cols_per_row)
                            for col, item in zip(row_cols, row_items):
                                with col:
                                    is_major = item["index"] == major_index
                                    with st.container(border=is_major):
                                        label = f"**Microspecies {item['index']}**"
                                        if is_major:
                                            label += " ⭐ Major"
                                        st.markdown(label)
                                        micro_img = smiles_to_2d_image(item["smiles"])
                                        if micro_img:
                                            st.image(micro_img, width="stretch")
                                        else:
                                            st.caption("Could not render structure.")
                                        micro_sdf_path = smiles_to_3d_file(item["smiles"], f"microspecies_{item['index']}.sdf")
                                        if micro_sdf_path:
                                            st.download_button(
                                                "Download .sdf",
                                                read_file_content(micro_sdf_path),
                                                f"microspecies_{item['index']}.sdf",
                                                "chemical/x-mdl-sdfile",
                                                key=f"dl_micro_{item['index']}",
                                            )

                        # --- Microspecies distribution swept across the full pH range ---
                        st.subheader("Microspecies distribution across the pH range")
                        st.caption(
                            "The same microspecies breakdown as above, now swept across pH 0–14 "
                            "in steps of 0.5. Each microspecies gets its own color in the stacked "
                            "area chart; at every pH the stack still adds up to 100%. The dashed "
                            "line marks the selected target pH."
                        )
                        sweep_df = build_macrospecies_distribution_sweep(states, ph_min=0.0, ph_max=14.0, step=0.5)
                        st.altair_chart(plot_macrospecies_sweep(sweep_df, target_ph), use_container_width=True)
                        st.download_button(
                            label="Download pH sweep data (.csv)",
                            data=sweep_df.to_csv(index=False),
                            file_name="microspecies_ph_sweep.csv",
                            mime="text/csv",
                        )

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
                            structure_path = smiles_to_3d_file(smi, f"protonated_{i}.sdf")
                            img = smiles_to_2d_image(smi)

                            if structure_path and img:
                                c3, c4 = st.columns(2)
                                with c3:
                                    st.image(img, width="stretch")
                                with c4:
                                    st_molstar(structure_path, key=f"molstar_res_{i}")
                                    content = read_file_content(structure_path)
                                    st.download_button("Download", content, f"structure_{i}.sdf", "chemical/x-mdl-sdfile", key=f"dl_{i}")
                else:
                    # Single result
                    smi, label = results[0]
                    st.subheader(label)
                    st.code(smi, language="smiles")
                    structure_path = smiles_to_3d_file(smi, f"protonated.sdf")
                    img = smiles_to_2d_image(smi)

                    if structure_path and img:
                        c3, c4 = st.columns(2)
                        with c3:
                            st.image(img, width="stretch")
                        with c4:
                            st_molstar(structure_path, key=f"molstar_res_single")
                            content = read_file_content(structure_path)
                            st.download_button("Download", content, "protonated_structure.sdf", "chemical/x-mdl-sdfile")
