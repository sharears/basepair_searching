import streamlit as st
import pandas as pd
import gdown
import os
import numpy as np
from collections import defaultdict
import pyarrow.parquet as pq
from PIL import Image

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


# ==================================================
# STEP 1 — Page setup (ALWAYS first)
# ==================================================
st.set_page_config(
    page_title="RNA Base Pair Hydrogen Bond Explorer",
    layout="wide"
)

#st.title("RNA Base Pair Hydrogen Bond Explorer")
img = Image.open("assets/weird_bps.png")
st.image(img, use_container_width=True)
st.markdown(
    """
    <div style="
        font-family: Arial, sans-serif;
        font-size: 36px;
        color: orange;
        text-align: center;
        margin-top: 6px;
    ">
        RNA Base Pair Hydrogen Bond Explorer
    </div>
    """,
    unsafe_allow_html=True
)


# ==================================================
# STEP 2 — Initialize session state (MUST come early)
# ==================================================
if "data_loaded" not in st.session_state:
    st.session_state.data_loaded = False

if "df_bp" not in st.session_state:
    st.session_state.df_bp = None


# ==================================================
# Data loader (Parquet from Google Drive)
# ==================================================
@st.cache_data(show_spinner=True)
def load_data_from_gdrive():
    url = "https://drive.google.com/file/d/1P1G0HMJnIOmg2ElosOEXtAlsPPGL8Rf0/view?usp=sharing"

    file_id = url.split("/d/")[1].split("/")[0]
    download_url = f"https://drive.google.com/uc?id={file_id}"

    parquet_file = "data.parquet"
    gdown.download(download_url, parquet_file, quiet=True)

    df = pd.read_parquet(parquet_file)
    os.remove(parquet_file)

    # 🔥 RECREATE combined_hbond_* columns (CRITICAL)
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
# STEP 4 — Guard clause (CRITICAL)
# ==================================================
if not st.session_state.data_loaded:
    st.info("Click **Load database** to start.")
    st.stop()


# ==================================================
# STEP 5 — Safe zone: use the data
# ==================================================
df_bp = st.session_state.df_bp

st.success(f"Database loaded: {len(df_bp):,} base pairs")


# ==================================================
# Everything below here is normal app logic
# ==================================================

# ==================================================
# User inputs
# ==================================================

st.subheader("Search criteria")

bp = st.selectbox(
    "Select base pair",
    sorted(df_bp["base_pair"].unique()),
    key="bp_select"
)

hbonds_input = st.text_input(
    "Hydrogen bonds (comma-separated)",
    placeholder="e.g. O6-N3, N2-O2",
    key="hbonds_input"
)

hbonds = [h.strip() for h in hbonds_input.split(",") if h.strip()]

#if st.button("Search"):
#    st.write("Search logic goes here")


# ==================================================
# Run search
# ==================================================
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
