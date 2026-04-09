import streamlit as st
import pandas as pd
import gdown
import os
import re
import json
import numpy as np
from collections import defaultdict
import pyarrow.parquet as pq
from PIL import Image
import base64
import py3Dmol
from groq import Groq

# ==================================================
# LLM CONFIG
# ==================================================
MODEL_NAME = "llama-3.3-70b-versatile"


def get_groq_client():
    """
    Create Groq client using Streamlit secrets first, then environment variable.
    """
    api_key = None

    if "GROQ_API_KEY" in st.secrets:
        api_key = st.secrets["GROQ_API_KEY"]
    else:
        api_key = os.environ.get("GROQ_API_KEY")

    if not api_key:
        return None

    return Groq(api_key=api_key)


# ==================================================
# Global CSS (labels + selectbox text)
# ==================================================
st.markdown(
    """
    <style>
    [data-testid="stWidgetLabel"] * {
        font-size: 24px !important;
        font-family: Arial, sans-serif !important;
        font-weight: 600 !important;
    }

    div[role="combobox"] {
        font-size: 22px !important;
        font-family: Arial, sans-serif !important;
    }

    div[role="listbox"] * {
        font-size: 22px !important;
        font-family: Arial, sans-serif !important;
        line-height: 1.4 !important;
    }
    </style>
    """,
    unsafe_allow_html=True
)


# ==================================================
# Allowed atoms dictionary
# ==================================================
ALLOWED_ATOMS = {
    "A": ["N1", "N3", "N6", "N7", "C2", "C8", "O2'"],
    "G": ["N1", "N2", "N3", "O6", "N7", "C8", "O2'"],
    "C": ["N3", "N4", "O2", "C5", "C6", "O2'"],
    "U": ["N3", "O2", "O4", "C5", "C6", "O2'"],
}


# ==================================================
# LLM parsing helpers
# ==================================================
def normalize_base_pair(bp):
    """
    Normalize base pair string to canonical format like G-C
    """
    bp = bp.strip().upper().replace("/", "-")
    parts = bp.split("-")
    if len(parts) != 2:
        return bp
    return f"{parts[0]}-{parts[1]}"


def parse_hbond_description_with_llm(bp, user_text):
    """
    Convert user natural-language description into strict hydrogen-bond list.
    Returns a list like:
    ["G.O6-C.N4", "G.N1-C.N3", "G.N2-C.O2"]
    """
    client = get_groq_client()
    if client is None:
        return []

    bp = normalize_base_pair(bp)

    prompt = f"""
You are helping parse RNA base-pair hydrogen bond requests.

The selected base pair is: {bp}

Allowed atoms:
A: {ALLOWED_ATOMS["A"]}
G: {ALLOWED_ATOMS["G"]}
C: {ALLOWED_ATOMS["C"]}
U: {ALLOWED_ATOMS["U"]}

Return ONLY valid JSON in this exact format:
{{
  "hbonds": ["X.ATOM-Y.ATOM", "X.ATOM-Y.ATOM"]
}}

Rules:
1. Use exactly the format Residue.Atom-Residue.Atom
2. Example: G.O6-C.N4
3. Do not include explanations
4. Do not include markdown
5. Do not include any text outside JSON
6. If the request is ambiguous or impossible, return:
   {{"hbonds": []}}
7. Use residue identities that are consistent with the selected base pair: {bp}

User request:
\"\"\"{user_text}\"\"\"
"""

    try:
        response = client.chat.completions.create(
            model=MODEL_NAME,
            temperature=0,
            messages=[
                {
                    "role": "system",
                    "content": "You convert RNA hydrogen-bond descriptions into strict JSON."
                },
                {
                    "role": "user",
                    "content": prompt
                }
            ]
        )

        content = response.choices[0].message.content.strip()
        parsed = json.loads(content)
        hbonds = parsed.get("hbonds", [])

        if not isinstance(hbonds, list):
            return []

        return [str(x).strip() for x in hbonds if str(x).strip()]

    except Exception:
        return []


def validate_hbonds(hbonds, bp):
    """
    Validate and clean LLM-generated or manually entered hydrogen bonds.
    Expected format: G.O6-C.N4
    """
    bp = normalize_base_pair(bp)
    bp_parts = bp.split("-")
    if len(bp_parts) != 2:
        return []

    allowed_residues = set(bp_parts)
    pattern = re.compile(r"^[ACGU]\.[A-Za-z0-9']+-[ACGU]\.[A-Za-z0-9']+$")

    cleaned = []

    for h in hbonds:
        h = h.strip()
        if not pattern.match(h):
            continue

        left, right = h.split("-")
        res1, atom1 = left.split(".", 1)
        res2, atom2 = right.split(".", 1)

        if res1 not in allowed_residues or res2 not in allowed_residues:
            continue

        if atom1 not in ALLOWED_ATOMS.get(res1, []):
            continue

        if atom2 not in ALLOWED_ATOMS.get(res2, []):
            continue

        cleaned.append(h)

    # Remove duplicates while preserving order
    seen = set()
    unique_cleaned = []
    for h in cleaned:
        if h not in seen:
            unique_cleaned.append(h)
            seen.add(h)

    return unique_cleaned


def parse_manual_hbonds(bp, manual_text):
    """
    Parse comma-separated manual input:
    G.O6-C.N4, G.N1-C.N3, G.N2-C.O2
    """
    raw = [x.strip() for x in manual_text.split(",") if x.strip()]
    return validate_hbonds(raw, bp)


# ==================================================
# Efficient hydrogen-bond matching (vectorized)
# ==================================================
def has_hbond(df, hbond):
    variants = [hbond, "-".join(hbond.split("-")[::-1])]
    mask = pd.Series(False, index=df.index)

    for i in range(1, 11):
        col = f"combined_hbond_{i}"
        if col in df.columns:
            col_data = df[col].astype(str)
            local_mask = pd.Series(False, index=df.index)
            for v in variants:
                local_mask |= col_data.str.contains(re.escape(v), na=False)
            mask |= local_mask

    return mask


def find_bp_interest(df, bp, hbonds):
    """
    Search dataframe by base pair and exact hydrogen-bond requirements.
    """
    bp_split = bp.split('-')
    if len(bp_split) != 2:
        return pd.DataFrame()

    # Normalize base pair orientation (GU == UG)
    if bp_split[0] != bp_split[1]:
        bps = [bp, "-".join(bp_split[::-1])]
    else:
        bps = [bp]

    mask = df["base_pair"].isin(bps)

    for hbond in hbonds:
        mask &= has_hbond(df, hbond)

    results = df[mask].copy()

    def extract_bp(row, hbond):
        variants = [hbond, "-".join(hbond.split("-")[::-1])]
        for i in range(1, 11):
            col = f'combined_hbond_{i}'
            if col in row.index:
                val = row[col]
                if isinstance(val, str):
                    for v in variants:
                        if v in val:
                            try:
                                return float(val.split('_')[-1])
                            except Exception:
                                return None
        return None

    for hbond in hbonds:
        results[hbond] = results.apply(lambda row: extract_bp(row, hbond), axis=1)

    col_drop = [col for col in results.columns if col.startswith('combined_hbond_')]
    col_drop += [col for col in results.columns if col.startswith('distance_hbond_')]
    col_drop += [col for col in results.columns if col.startswith('atoms_ID_hbond_')]

    drop_cols = ['nt1', 'nt2'] + col_drop
    drop_cols = [c for c in drop_cols if c in results.columns]

    results1 = results.drop(columns=drop_cols)
    return results1


# ==================================================
# 3D structure rendering helper
# ==================================================
def render_basepair_3d(
    pdb_id,
    chain1, resi1, icode1,
    chain2, resi2, icode2,
    base1=None, base2=None
):
    view = py3Dmol.view(query=f"pdb:{str(pdb_id).lower()}", width=700, height=520)
    view.setBackgroundColor("white")

    base_colors = {"A": "orange", "G": "green", "C": "cyan", "U": "brown"}

    def norm_chain(x):
        return str(x).strip()

    def norm_resi(x):
        return int(float(x))

    def norm_icode(x):
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return None
        s = str(x).strip()
        if not s or s.lower() == "nan":
            return None
        return s

    def make_sel(chain, resi, icode):
        sel = {"chain": norm_chain(chain), "resi": norm_resi(resi)}
        ic = norm_icode(icode)
        if ic:
            sel["icode"] = ic
        return sel

    sel1 = make_sel(chain1, resi1, icode1)
    sel2 = make_sel(chain2, resi2, icode2)

    view.setStyle({}, {"stick": {"radius": 0.12, "color": "lightgray"}})
    view.zoomTo()

    view.setStyle(
        sel1,
        {"stick": {"radius": 0.30, "color": base_colors.get(str(base1).strip(), "orange")}}
    )
    view.setStyle(
        sel2,
        {"stick": {"radius": 0.30, "color": base_colors.get(str(base2).strip(), "green")}}
    )

    view.zoomTo({"or": [sel1, sel2]})

    label1 = f"{norm_chain(chain1)}:{norm_resi(resi1)}{norm_icode(icode1) or ''}"
    label2 = f"{norm_chain(chain2)}:{norm_resi(resi2)}{norm_icode(icode2) or ''}"

    view.addLabel(
        f"{str(pdb_id).upper()}  {label1} – {label2}",
        {
            "fontColor": "black",
            "backgroundColor": "white",
            "borderColor": "gray",
            "borderThickness": 1
        }
    )

    st.components.v1.html(view._make_html(), height=540)


# ==================================================
# STEP 1 — Page setup
# ==================================================
def svg_to_base64(svg_path):
    with open(svg_path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")


st.set_page_config(layout="wide")

svg_b64 = svg_to_base64("assets/weird_bps.svg")

st.markdown(
    f"""
    <div style="display:flex; justify-content:center;">
        <img src="data:image/svg+xml;base64,{svg_b64}"
             style="width:70vw; max-width:1400px;">
    </div>
    """,
    unsafe_allow_html=True
)

st.markdown(
    """
    <div style="
        font-family: Arial, sans-serif;
        font-size: 36px;
        color: orange;
        text-align: center;
        margin-top: 12px;
    ">
        RNA Base Pair Explorer
    </div>
    """,
    unsafe_allow_html=True
)


# ==================================================
# STEP 2 — Initialize session state
# ==================================================
if "data_loaded" not in st.session_state:
    st.session_state.data_loaded = False

if "df_bp" not in st.session_state:
    st.session_state.df_bp = None

if "results" not in st.session_state:
    st.session_state.results = None

if "parsed_hbonds" not in st.session_state:
    st.session_state.parsed_hbonds = []

if "search_error" not in st.session_state:
    st.session_state.search_error = None


# ==================================================
# Data loader (Parquet from Google Drive)
# ==================================================
@st.cache_data(show_spinner=True)
def load_data_from_gdrive():
    url = 'https://drive.google.com/file/d/1MH0Xjo8RjJMdIoUbpQlUEWrZfpMeG2yX/view?usp=sharing'

    file_id = url.split("/d/")[1].split("/")[0]
    download_url = f"https://drive.google.com/uc?id={file_id}"

    parquet_file = "data.parquet"
    gdown.download(download_url, parquet_file, quiet=True)

    df = pd.read_parquet(parquet_file)
    os.remove(parquet_file)

    suffix_groups = defaultdict(list)
    for col in df.columns:
        suffix = "_".join(col.split("_")[-2:])
        if suffix.startswith("hbond"):
            suffix_groups[suffix].append(col)

    for suffix, cols in suffix_groups.items():
        if len(cols) == 2:
            atom_col, dist_col = sorted(cols)
            df[f"combined_{suffix}"] = (
                df[atom_col].astype(str) + "_" +
                df[dist_col].astype(str)
            )

    return df


# ==================================================
# STEP 3 — User-triggered heavy action
# ==================================================
if st.button("Load database"):
    with st.spinner("Loading database..."):
        st.session_state.df_bp = load_data_from_gdrive()
        st.session_state.data_loaded = True


# ==================================================
# STEP 4 — Guard clause
# ==================================================
if not st.session_state.data_loaded:
    st.info("Click **Load database** to start.")
    st.stop()


# ==================================================
# STEP 5 — Safe zone
# ==================================================
df_bp = st.session_state.df_bp

st.markdown(
    f"""
    <div style="
        background-color: #e6f4ea;
        border-left: 6px solid #2e7d32;
        padding: 14px 18px;
        border-radius: 6px;
        font-family: Arial, sans-serif;
        font-size: 22px;
        color: #1b5e20;
        margin-top: 10px;
        margin-bottom: 10px;
    ">
        <b>Database loaded:</b> {df_bp['PDB_ID'].nunique():,} unique structures
        ({len(df_bp):,} base pairs)
    </div>
    """,
    unsafe_allow_html=True
)

st.markdown(
    """
    <style>
    span[data-testid="stWidgetLabel"] {
        font-size: 30px !important;
        font-family: Arial, sans-serif !important;
        font-weight: 600;
    }
    </style>
    """,
    unsafe_allow_html=True
)


# ==================================================
# User inputs
# ==================================================
st.markdown(
    """
    <div style="
        font-family: Arial, sans-serif;
        font-size: 30px;
        font-weight: 700;
        color: white;
        margin-bottom: 12px;
    ">
        Search criteria
    </div>
    """,
    unsafe_allow_html=True
)

col1, col2 = st.columns([1, 2])

with col1:
    bp = st.selectbox(
        "Select base pair",
        sorted(df_bp["base_pair"].unique()),
        key="bp_select"
    )

with col2:
    input_mode = st.radio(
        "Hydrogen-bond input mode",
        ["Natural language", "Manual list"],
        horizontal=True,
        key="input_mode"
    )

    if input_mode == "Natural language":
        hbond_description = st.text_area(
            "Describe the hydrogen bonds you want",
            placeholder=(
                "Example: standard Watson–Crick G-C hydrogen bonds, "
                "with G.O6 to C.N4, G.N1 to C.N3, and G.N2 to C.O2"
            ),
            key="hbond_description",
            height=120
        )
        manual_hbonds_input = ""
    else:
        manual_hbonds_input = st.text_input(
            "Hydrogen bonds (comma-separated)",
            placeholder="e.g. G.O6-C.N4, G.N1-C.N3, G.N2-C.O2",
            key="hbonds_input"
        )
        hbond_description = ""


# ==================================================
# Run search
# ==================================================
if st.button("Search"):
    st.session_state.search_error = None
    st.session_state.parsed_hbonds = []
    st.session_state.results = None

    if input_mode == "Natural language":
        if not hbond_description.strip():
            st.session_state.search_error = "Please describe at least one hydrogen bond."
        else:
            with st.spinner("Interpreting hydrogen-bond description..."):
                parsed_hbonds = parse_hbond_description_with_llm(bp, hbond_description)
                parsed_hbonds = validate_hbonds(parsed_hbonds, bp)

            if not parsed_hbonds:
                st.session_state.search_error = (
                    "Could not interpret a valid hydrogen-bond pattern from your description."
                )
            else:
                st.session_state.parsed_hbonds = parsed_hbonds
                with st.spinner("Searching base-pair database..."):
                    st.session_state.results = find_bp_interest(df_bp, bp, parsed_hbonds)

    else:
        parsed_hbonds = parse_manual_hbonds(bp, manual_hbonds_input)

        if not parsed_hbonds:
            st.session_state.search_error = (
                "Please enter valid hydrogen bonds like G.O6-C.N4, G.N1-C.N3."
            )
        else:
            st.session_state.parsed_hbonds = parsed_hbonds
            with st.spinner("Searching base-pair database..."):
                st.session_state.results = find_bp_interest(df_bp, bp, parsed_hbonds)


# ==================================================
# Display parse/search status
# ==================================================
if st.session_state.search_error:
    st.warning(st.session_state.search_error)

if st.session_state.parsed_hbonds:
    st.markdown(
        f"**Interpreted hydrogen bonds:** {', '.join(st.session_state.parsed_hbonds)}"
    )


# ==================================================
# Results
# ==================================================
results = st.session_state.results

if results is not None:
    st.markdown(f"### Results ({len(results)} matches)")

    if results.empty:
        st.info("No matching base pairs found.")
    else:
        st.dataframe(results.head(1000), use_container_width=True)

        # ==================================================
        # Interactive 3D base-pair viewer
        # ==================================================
        st.subheader("Interactive 3D base-pair view")

        selected_idx = st.selectbox(
            "Select base-pair instance to visualize",
            results.index,
            format_func=lambda i: (
                f"{results.loc[i, 'PDB_ID']} | "
                f"{results.loc[i, 'chain_ID_res1']}:{results.loc[i, 'res_index_res1']} "
                f"{results.loc[i, 'res_ID_res1']} — "
                f"{results.loc[i, 'chain_ID_res2']}:{results.loc[i, 'res_index_res2']} "
                f"{results.loc[i, 'res_ID_res2']}"
            ),
            key="bp_3d_select"
        )

        if st.button("Show 3D structure", key="show_3d"):
            row = results.loc[selected_idx]

            st.caption(
                f"3D selections: "
                f"sel1={row['chain_ID_res1']}:{row['res_index_res1']} icode={row.get('icode_res1')} | "
                f"sel2={row['chain_ID_res2']}:{row['res_index_res2']} icode={row.get('icode_res2')}"
            )

            render_basepair_3d(
                pdb_id=row["PDB_ID"],
                chain1=row["chain_ID_res1"],
                resi1=row["res_index_res1"],
                icode1=row.get("icode_res1"),
                chain2=row["chain_ID_res2"],
                resi2=row["res_index_res2"],
                icode2=row.get("icode_res2"),
                base1=row["res_ID_res1"],
                base2=row["res_ID_res2"]
            )

        st.caption(f"Showing first 1,000 of {len(results)} matches")

        st.download_button(
            label="Download full results as CSV",
            data=results.to_csv(index=False),
            file_name="bp_search_results.csv",
            mime="text/csv"
        )