import streamlit as st
import pandas as pd
import gdown
import os

def find_bp_interest(df, bp, hbonds):
    # Normalize base pair (GU == UG)
    bp_split = bp.split('-')
    if bp_split[0] != bp_split[1]:
        bps = [bp, "-".join(bp_split[::-1])]
    else:
        bps = [bp]

    df1 = df[df['base_pair'].isin(bps)].copy()
    df1.index = np.arange(len(df1))

    if df1.empty:
        return df1

    # Group hbond columns by suffix (e.g., hbond_1)
    suffix_groups = defaultdict(list)
    for col in df1.columns:
        suffix = '_'.join(col.split('_')[-2:])
        if suffix.startswith('hbond'):
            suffix_groups[suffix].append(col)

    # Combine atom + distance into combined_hbond_i
    for suffix, cols in suffix_groups.items():
        if len(cols) == 2:
            atom_col, dist_col = sorted(cols)
            new_col = f'combined_{suffix}'
            df1[new_col] = (
                df1[atom_col].astype(str) + "_" +
                df1[dist_col].astype(str)
            )

    # Helper: check if a hydrogen bond exists in a row
    def extract_bp(row, hbond):
        variants = [hbond, "-".join(hbond.split("-")[::-1])]

        for i in range(1, 11):  # up to 10 hbonds
            col = f'combined_hbond_{i}'
            if col in row:
                val = row[col]
                if isinstance(val, str):
                    for v in variants:
                        if v in val:
                            try:
                                return float(val.split('_')[-1])
                            except ValueError:
                                return None
        return None

    # Apply hydrogen bond filters
    for hbond in hbonds:
        df1[hbond] = df1.apply(
            lambda row: extract_bp(row, hbond), axis=1
        )

    # Keep rows containing ALL requested hbonds
    df2 = df1.dropna(subset=hbonds)

    # Clean up temporary columns
    df3 = df2.drop(
        columns=[c for c in df2.columns if c.startswith("combined_")],
        errors="ignore"
    )
    df3.index = np.arange(len(df3))

    return df3

# --------------------------------------------------
# Page setup
# --------------------------------------------------
st.set_page_config(
    page_title="RNA Base Pair Hydrogen Bond Explorer",
    layout="wide"
)

st.title("RNA Base Pair Hydrogen Bond Explorer")

# --------------------------------------------------
# Load data from Google Drive (cached)
# --------------------------------------------------
@st.cache_data(show_spinner=True)
def load_data_from_gdrive():
    #url = "https://drive.google.com/file/d/1aNb12ww1SF3ydfr3vqrdZ7yvpNK82voG/view?usp=drive_link"

    url= 'https://drive.google.com/file/d/1aNb12ww1SF3ydfr3vqrdZ7yvpNK82voG/view?usp=drive_link'

    # Extract file ID
    file_id = url.split('/d/')[1].split('/')[0]
    download_url = f"https://drive.google.com/uc?id={file_id}"

    output = "data.csv"

    # Download
    gdown.download(download_url, output, quiet=True)

    # Load CSV
    df = pd.read_csv(output)

    # Clean up
    os.remove(output)

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
            st.dataframe(results, use_container_width=True)

            st.download_button(
                label="Download results as CSV",
                data=results.to_csv(index=False),
                file_name="bp_search_results.csv",
                mime="text/csv"
            )
