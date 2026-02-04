import streamlit as st
import pandas as pd
import gdown
import os

from bp_search import find_bp_interest   # or paste function directly here

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
    url = "https://drive.google.com/file/d/1aNb12ww1SF3ydfr3vqrdZ7yvpNK82voG/view?usp=drive_link"

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
