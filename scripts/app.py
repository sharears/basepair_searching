import streamlit as st
import pandas as pd
import gdown
import os
import numpy as np
from collections import defaultdict

# --------------------------------------------------
# Efficient hydrogen-bond matching (vectorized)
# --------------------------------------------------
def has_hbond(df, hbond):
    variants = [hbond, "-".join(hbond.split("-")[::-1])]
    mask = False

    for i in range(1, 11):
        col = f"combined_hbond_{i}"
        if col in df.columns:
            for v in variants:
                mask |= df[col].str.contains(v, na=False)

    return mask


def find_bp_interest(df, bp, hbonds):
    # Normalize base pair (GU == UG)
    bp_split = bp.split('-')
    if bp_split[0] != bp_split[1]:
        bps = [bp, "-".join(bp_split[::-1])]
    else:
        bps = [bp]

    # 🔥 No copy here (saves memory)
    mask = df["base_pair"].isin(bps)

    # Apply hydrogen bond filters
    for hbond in hbonds:
        mask &= has_hbond(df, hbond)

    results = df[mask]

    return results


# --------------------------------------------------
# Page setup
# --------------------------------------------------
st.set_page_config(
    page_title="RNA Base Pair Hydrogen Bond Explorer",
    layout="wide"
)

st.title("RNA Base Pair Hydrogen Bond Explorer")

# --------------------------------------------------
# Load & preprocess data ONCE (cached)
# --------------------------------------------------
@st.cache_data(show_spinner=True)
def load_data_from_gdrive():
    url = "https://drive.google.com/file/d/1aNb12ww1SF3ydfr3vqrdZ7yvpNK82voG/view?usp=drive_link"

    # Sample file for quick testing
    # url = "https://drive.google.com/file/d/1NWrU7Xu-7K-5wRE1gdIpUIFNe_6jhrRI/view?usp=drive_link"

    file_id = url.split('/d/')[1].split('/')[0]
    download_url = f"https://drive.google.com/uc?id={file_id}"

    output = "data.csv"
    gdown.download(download_url, output, quiet=True)

    # Load CSV (avoid dtype warning + extra memory churn)
    df = pd.read_csv(output, low_memory=False)
    os.remove(output)

    # 🔥 PRECOMPUTE combined_hbond columns ONCE
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


df_bp = load_data_from_gdrive()
st.success(f"Database loaded: {len(df_bp):,} base pairs")

# --------------------------------------------------
# User inputs
# --------------------------------------------------
st.subheader("Search criteria")

col1, col2 = st.columns(2)

with col1:
    bp = st.selectbox(
        "Select base pair",
        sorted(df_bp["base_pair"].unique())
    )

with col2:
    hbonds_input = st.text_input(
        "Hydrogen bonds (comma-separated)",
        placeholder="e.g. O6-N3, N2-O2"
    )

hbonds = [h.strip() for h in hbonds_input.split(",") if h.strip()]

# --------------------------------------------------
# Run search
# --------------------------------------------------
if st.button("Search"):
    if not hbonds:
        st.warning("Please enter at least one hydrogen bond.")
    else:
        with st.spinner("Searching base-pair database..."):
            results = find_bp_interest(df_bp, bp, hbonds)

        st.markdown(f"### Results ({len(results)} matches)")

        if results.empty:
            st.info("No matching base pairs found.")
        else:
            # 🔥 Do NOT render huge tables
            st.dataframe(results.head(1000), use_container_width=True)
            st.caption(
                f"Showing first 1,000 of {len(results)} matches"
            )

            st.download_button(
                label="Download full results as CSV",
                data=results.to_csv(index=False),
                file_name="bp_search_results.csv",
                mime="text/csv"
            )
