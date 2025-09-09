import streamlit as st
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
from io import StringIO
#
# --- Sequence Functions ---
def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return round(((g + c) / len(seq)) * 100, 2)

def codon_usage(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    return Counter(codons)

def find_orfs(seq, min_len=100):
    seq = seq.upper()
    is_rna = "U" in seq and "T" not in seq

    if is_rna:
        start_codon = "AUG"
        stop_codons = ["UAA", "UAG", "UGA"]
    else:
        start_codon = "ATG"
        stop_codons = ["TAA", "TAG", "TGA"]

    orfs = []
    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            codon = seq[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(seq)-2, 3):
                    if seq[j:j+3] in stop_codons:
                        length = j+3 - i
                        if length >= min_len:
                            orfs.append({
                                "frame": frame+1,
                                "start": i,
                                "end": j+3,
                                "length": length,
                                "type": "RNA" if is_rna else "DNA"
                            })
                        break
            i += 3
    return orfs

# --- Streamlit App ---
st.title("DNA/RNA Sequence Analyzer")
st.markdown(
    """
**DNA/RNA Sequence Analyzer** – An interactive Streamlit app in Python that automatically detects and analyzes DNA/RNA sequences from FASTA files.  
Features include GC content calculation, codon usage profiling with visualizations, and open reading frame (ORF) detection.  
Deployed via GitHub and Streamlit Cloud.
"""
)

use_example = st.checkbox("Use example sequence")

if uploaded or use_example:
    if use_example:
        seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAA"
        record_id = "Example_Sequence"
    else:
        stringio = StringIO(uploaded.getvalue().decode("utf-8"))
        record = next(SeqIO.parse(stringio, "fasta"))
        seq = str(record.seq).upper()
        record_id = record.id

    st.subheader("Sequence Info")
    st.write(f"ID: {record_id}")
    st.write(f"Length: {len(seq)} bp")
    st.write(f"GC Content: {gc_content(seq)} %")


uploaded = st.file_uploader("Upload a FASTA file", type=["fasta", "fa"])
if uploaded:
    # Convert uploaded file into StringIO for SeqIO
    stringio = StringIO(uploaded.getvalue().decode("utf-8"))
    record = next(SeqIO.parse(stringio, "fasta"))
    seq = str(record.seq).upper()
    
    st.subheader("Sequence Info")
    st.write(f"ID: {record.id}")
    st.write(f"Length: {len(seq)} bp")
    st.write(f"GC Content: {gc_content(seq)} %")
    
    st.subheader("Codon Usage")
    codons = codon_usage(seq)
    codon_df = {k: v for k, v in codons.items() if len(k) == 3}
    
    fig, ax = plt.subplots(figsize=(10,4))
    ax.bar(codon_df.keys(), codon_df.values())
    ax.set_title("Codon Usage")
    ax.set_xlabel("Codon")
    ax.set_ylabel("Frequency")
    plt.xticks(rotation=90)
    st.pyplot(fig)

    st.subheader("Open Reading Frames (≥100 bp)")
    orfs = find_orfs(seq)
    if not orfs:
        st.write("No ORFs found ≥100 bp.")
    else:
        for orf in orfs:
            st.write(
                f"Frame {orf['frame']}: {orf['start']}-{orf['end']} "
                f"({orf['length']} bp, {orf['type']})"
            )
