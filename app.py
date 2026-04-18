#!/usr/bin/env python3
"""
MetalCytoToxDB — Local Mirror
ML-assisted cytotoxicity database for transition metal complexes.
Based on Zenodo record 15853577.
"""

import streamlit as st
import sqlite3
import os
import io
import csv
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

# Setup paths
APP_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(APP_DIR, "data")
DB_PATH = os.path.join(DATA_DIR, "biometaldb.sqlite")
MOL_DIR = os.path.join(DATA_DIR, "mol")
MOL3_DIR = os.path.join(DATA_DIR, "mol3")
FAVICON_PATH = os.path.join(APP_DIR, "favicon_64.png")

# Add patches to path
sys.path.insert(0, os.path.join(APP_DIR, "patches"))
from selfies_metal import selfies_to_smiles_metal

# Page config
st.set_page_config(
    page_title="MetalCytoToxDB",
    page_icon=FAVICON_PATH,
    layout="wide",
    initial_sidebar_state="expanded",
)

# Sidebar
st.sidebar.markdown("## DATABASE")
st.sidebar.markdown("### MetalCytoToxDB")
st.sidebar.markdown("---")

page = st.sidebar.radio("Navigation", [
    "🔍 Search complexes",
    "📚 Literature",
    "☀️ Phototoxicity",
    "📊 Statistics",
    "⚗️ Structures",
])

st.sidebar.markdown("---")
st.sidebar.markdown("""
**Cite this work:**  
Krasnov et al.  
*Machine Learning for Anticancer Activity Prediction of Transition Metal Complexes*  
[ChemRxiv 2025](https://doi.org/10.26434/chemrxiv-2025-1nqvm-v2) ↗

[Zenodo Dataset ↗](https://doi.org/10.5281/zenodo.15853577)
""")


@st.cache_resource
def get_db():
    """Get database connection."""
    return sqlite3.connect(DB_PATH, check_same_thread=False)


@st.cache_data
def get_stats():
    """Get database statistics."""
    conn = get_db()
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM complexes")
    n_complexes = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM measurements")
    n_measurements = cur.fetchone()[0]
    cur.execute("SELECT COUNT(DISTINCT cell_line) FROM measurements")
    n_cells = cur.fetchone()[0]
    cur.execute("SELECT COUNT(DISTINCT doi) FROM measurements")
    n_sources = cur.fetchone()[0]
    cur.execute("SELECT metal, COUNT(*) FROM complexes GROUP BY metal ORDER BY COUNT(*) DESC")
    metal_dist = cur.fetchall()
    return n_complexes, n_measurements, n_cells, n_sources, metal_dist


def render_mol_image(smiles, size=(400, 300)):
    """Render SMILES as 2D image."""
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=size)
        return img
    except Exception:
        return None


# ==================== SEARCH PAGE ====================
if page == "🔍 Search complexes":
    st.title("ML-ASSISTED CYTOTOXICITY DATABASE")
    st.header("MetalCytoToxDB")
    st.markdown("Explore IC₅₀ cytotoxicity data for transition metal complexes. Search by cell line, ligand structure, or literature source.")
    
    # Stats header
    n_complexes, n_measurements, n_cells, n_sources, metal_dist = get_stats()
    
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("IC₅₀ VALUES", f"{n_measurements:,}")
    col2.metric("COMPLEXES", f"{n_complexes:,}")
    col3.metric("CELL LINES", f"{n_cells:,}")
    col4.metric("SOURCES", f"{n_sources:,}")
    
    # Metal distribution
    st.markdown("---")
    cols = st.columns(len(metal_dist))
    for i, (metal, count) in enumerate(metal_dist):
        pct = count / n_complexes * 100
        cols[i].metric(metal, f"{count:,}", f"{pct:.1f}%")
    
    # Search
    st.markdown("---")
    smiles_input = st.text_input("SMILES", placeholder="Enter SMILES string for ligand search...")
    
    # Popular ligands
    st.markdown("**Popular ligands:**")
    ligand_cols = st.columns(7)
    popular = ["p-cymene", "ppy", "bpy", "PPh3", "phen", "Cp", "dfppy"]
    for i, lig in enumerate(popular):
        if ligand_cols[i].button(lig):
            smiles_input = lig
    
    # Search mode
    search_mode = st.radio("Search mode", [
        "Full molecule match", "Substructure search", "Similarity search"
    ], horizontal=True)
    
    # Filters
    conn = get_db()
    cur = conn.cursor()
    
    filter_cols = st.columns(6)
    metal_filter = filter_cols[0].selectbox("METAL", ["All metals"] + [m for m, _ in metal_dist])
    
    # Donor atom filter
    cur.execute("""
        SELECT donor_atoms, COUNT(*) as cnt 
        FROM complexes 
        WHERE donor_atoms IS NOT NULL AND donor_atoms != 'None' AND donor_atoms != ''
        GROUP BY donor_atoms 
        ORDER BY cnt DESC 
        LIMIT 20
    """)
    donor_options = ["All donors"] + [f"{da} ({cnt})" for da, cnt in cur.fetchall()]
    donor_filter = filter_cols[1].selectbox("DONOR ATOMS", donor_options)
    
    # Cell lines
    cur.execute("SELECT DISTINCT cell_line FROM measurements ORDER BY cell_line")
    cell_lines = [r[0] for r in cur.fetchall() if r[0]]
    cell_filter = filter_cols[2].selectbox("CELL LINE", ["All cell lines"] + cell_lines)
    
    time_filter = filter_cols[3].selectbox("EXPOSURE TIME", ["All time ranges", "24h", "48h", "72h"])
    ic50_min = filter_cols[4].number_input("IC₅₀ min, μM", min_value=0.0, value=0.0, step=0.01)
    ic50_max = filter_cols[5].number_input("IC₅₀ max, μM", min_value=0.0, value=0.0, step=0.01)
    
    # Build query
    query = """
        SELECT c.id, c.metal, c.oxidation_state, c.smiles_ligands, c.metal_smiles, c.donor_atoms,
               m.cell_line, m.ic50_dark_value, m.exposure_time_h, m.doi, m.year, m.ic50_cisplatin_value
        FROM measurements m
        JOIN complexes c ON m.complex_id = c.id
        WHERE 1=1
    """
    params = []
    
    if metal_filter != "All metals":
        query += " AND c.metal = ?"
        params.append(metal_filter)
    
    if donor_filter != "All donors":
        donor_val = donor_filter.split(" (")[0]  # Extract JSON from "{'N': 1} (783)"
        query += " AND c.donor_atoms = ?"
        params.append(donor_val)
    
    if cell_filter != "All cell lines":
        query += " AND m.cell_line = ?"
        params.append(cell_filter)
    
    if smiles_input and smiles_input not in popular:
        if search_mode == "Full molecule match":
            query += " AND c.smiles_ligands = ?"
            params.append(smiles_input)
        elif search_mode == "Substructure search":
            # Simple contains search
            query += " AND c.smiles_ligands LIKE ?"
            params.append(f"%{smiles_input}%")
    
    if ic50_min > 0:
        query += " AND m.ic50_dark_value >= ?"
        params.append(ic50_min)
    if ic50_max > 0:
        query += " AND m.ic50_dark_value <= ?"
        params.append(ic50_max)
    
    if time_filter != "All time ranges":
        hours = float(time_filter.replace("h", ""))
        query += " AND m.exposure_time_h = ?"
        params.append(hours)
    
    # Sorting
    query += " ORDER BY m.ic50_dark_value ASC NULLS LAST LIMIT 500"
    
    cur.execute(query, params)
    results = cur.fetchall()
    
    st.markdown(f"### RESULTS — {len(results)} records")
    
    if results:
        # Group by complex
        from collections import defaultdict
        complex_groups = defaultdict(list)
        for r in results:
            cid = r[0]
            complex_groups[cid].append(r)
        
        st.markdown(f"**SHOWING** 1–{min(len(complex_groups), 20)} of {len(complex_groups)} complexes")
        
        # Download CSV
        csv_buf = io.StringIO()
        writer = csv.writer(csv_buf)
        writer.writerow(["Metal", "Oxidation", "SMILES", "Donor Atoms", "Cell Line", "IC50 (μM)", "Cisplatin IC50 (μM)", "DOI", "Year"])
        for r in results:
            writer.writerow([r[1], r[2], r[3], r[5], r[6], r[7], r[11], r[9], r[10]])
        
        st.download_button("Download CSV", csv_buf.getvalue(), "metalcytotoxdb_results.csv", "text/csv")
        
        # Display complexes
        for i, (cid, group) in enumerate(list(complex_groups.items())[:20]):
            r = group[0]
            metal = r[1]
            ox = r[2]
            smiles = r[3]
            metal_smiles = r[4]
            donor_atoms = r[5]
            
            min_ic50 = min((g[7] for g in group if g[7]), default=None)
            n_cells = len(set(g[6] for g in group if g[6]))
            
            donor_str = f" | donors: {donor_atoms}" if donor_atoms and donor_atoms != "None" else ""
            with st.expander(f"{metal}({ox}){donor_str} — {min_ic50:.2g} μM min IC₅₀ · {n_cells} cell lines" if min_ic50 else f"{metal}({ox}){donor_str} — {n_cells} cell lines"):
                # Render structure
                img = render_mol_image(metal_smiles or smiles)
                if img:
                    st.image(img, width=400)
                
                st.code(smiles, language="text")
                
                # MOL download
                mol_path = os.path.join(MOL_DIR, f"complex_{cid}.mol")
                if os.path.exists(mol_path):
                    with open(mol_path, 'r') as f:
                        st.download_button(f"Download MOL (V2000)", f.read(), f"complex_{cid}.mol", "chemical/x-mdl-molfile", key=f"mol_{cid}")
                
                mol3_path = os.path.join(MOL3_DIR, f"complex_{cid}.mol")
                if os.path.exists(mol3_path):
                    with open(mol3_path, 'r') as f:
                        st.download_button(f"Download MOL (V3000)", f.read(), f"complex_{cid}_v3000.mol", "chemical/x-mdl-molfile", key=f"mol3_{cid}")
                
                # IC50 table
                st.markdown("**CELL LINE** | **TIME** | **IC₅₀, μM** | **CISPLATIN** | **DOI**")
                for g in group[:10]:
                    cl = g[6] or ""
                    t = f"{g[8]:.0f}h" if g[8] else ""
                    ic50 = f"{g[7]:.2g}" if g[7] else "—"
                    cis = f"{g[11]:.2g}" if g[11] else "—"
                    doi = g[9] or ""
                    st.markdown(f"{cl} | {t} | {ic50} | {cis} | [{doi}](https://doi.org/{doi})")


# ==================== LITERATURE PAGE ====================
elif page == "📚 Literature":
    st.title("📚 Literature")
    
    conn = get_db()
    cur = conn.cursor()
    
    # Search/filter
    search_query = st.text_input("Search papers by title, author, or DOI", "")
    
    # Build query with paper metadata
    query = """
        SELECT p.doi, p.year, p.title, p.authors, 
               COUNT(DISTINCT m.complex_id) as n_complexes,
               COUNT(*) as n_records,
               GROUP_CONCAT(DISTINCT c.metal) as metals
        FROM measurements m
        JOIN complexes c ON m.complex_id = c.id
        LEFT JOIN papers p ON m.doi = p.doi
        WHERE m.doi IS NOT NULL AND m.doi != ''
    """
    params = []
    
    if search_query:
        query += " AND (p.title LIKE ? OR p.authors LIKE ? OR m.doi LIKE ?)"
        params.extend([f"%{search_query}%", f"%{search_query}%", f"%{search_query}%"])
    
    query += " GROUP BY m.doi ORDER BY n_records DESC"
    
    if not search_query:
        query += " LIMIT 200"
    
    cur.execute(query, params)
    papers = cur.fetchall()
    
    st.markdown(f"### {len(papers)} papers" + (f" matching '{search_query}'" if search_query else " (top by data points)"))
    
    # Table view
    if papers:
        for doi, year, title, authors, n_complexes, n_records, metals in papers:
            title_display = title if title else "(title pending)"
            authors_short = (authors[:80] + "...") if authors and len(authors) > 80 else (authors or "")
            year_str = str(year) if year else "?"
            
            with st.expander(f"{year_str} | {n_records} records | {title_display[:80]}"):
                st.markdown(f"**Title:** {title or 'N/A'}")
                st.markdown(f"**Authors:** {authors or 'N/A'}")
                st.markdown(f"**Year:** {year_str}")
                st.markdown(f"**Metals:** {metals or 'N/A'}")
                st.markdown(f"**Complexes:** {n_complexes} | **Measurements:** {n_records}")
                st.markdown(f"**DOI:** [{doi}](https://doi.org/{doi})")


# ==================== PHOTOTOXICITY PAGE ====================
elif page == "☀️ Phototoxicity":
    st.title("☀️ Phototoxicity")
    
    conn = get_db()
    cur = conn.cursor()
    
    cur.execute("""
        SELECT c.metal, c.smiles_ligands, m.cell_line, 
               m.ic50_dark_value, m.ic50_light_value, m.doi
        FROM measurements m
        JOIN complexes c ON m.complex_id = c.id
        WHERE m.ic50_dark_value IS NOT NULL AND m.ic50_light_value IS NOT NULL
        ORDER BY m.ic50_dark_value ASC
        LIMIT 200
    """)
    results = cur.fetchall()
    
    st.markdown(f"### {len(results)} complexes with both dark and light IC₅₀ data")
    
    if results:
        st.markdown("**Metal** | **Cell Line** | **Dark IC₅₀, μM** | **Light IC₅₀, μM** | **PI**")
        for metal, smiles, cl, dark, light, doi in results[:50]:
            pi = dark / light if light and light > 0 else "—"
            pi_str = f"{pi:.2f}" if isinstance(pi, float) else pi
            st.markdown(f"{metal} | {cl} | {dark:.2g} | {light:.2g} | {pi_str}")


# ==================== STATISTICS PAGE ====================
elif page == "📊 Statistics":
    st.title("📊 Statistics")
    
    n_complexes, n_measurements, n_cells, n_sources, metal_dist = get_stats()
    
    st.metric("Total IC₅₀ measurements", f"{n_measurements:,}")
    st.metric("Unique complexes", f"{n_complexes:,}")
    st.metric("Cell lines", f"{n_cells:,}")
    st.metric("Literature sources", f"{n_sources:,}")
    
    st.markdown("### Metal distribution")
    for metal, count in metal_dist:
        pct = count / n_complexes * 100
        st.progress(pct / 100, text=f"{metal}: {count:,} ({pct:.1f}%)")
    
    # Year distribution
    conn = get_db()
    cur = conn.cursor()
    cur.execute("SELECT year, COUNT(*) FROM measurements WHERE year IS NOT NULL GROUP BY year ORDER BY year")
    year_data = cur.fetchall()
    
    if year_data:
        st.markdown("### Publications by year")
        st.bar_chart({str(y): c for y, c in year_data})
    
    # Cell line distribution
    cur.execute("SELECT cell_line, COUNT(*) as n FROM measurements WHERE cell_line != '' GROUP BY cell_line ORDER BY n DESC LIMIT 20")
    cell_data = cur.fetchall()
    if cell_data:
        st.markdown("### Top 20 cell lines")
        st.bar_chart({cl: n for cl, n in cell_data})
    
    # Donor atom distribution
    cur.execute("""
        SELECT donor_atoms, COUNT(*) as cnt 
        FROM complexes 
        WHERE donor_atoms IS NOT NULL AND donor_atoms != 'None' AND donor_atoms != ''
        GROUP BY donor_atoms 
        ORDER BY cnt DESC 
        LIMIT 15
    """)
    donor_data = cur.fetchall()
    if donor_data:
        st.markdown("### Top 15 donor atom patterns")
        st.bar_chart({da: cnt for da, cnt in donor_data})
    
    # Oxidation state distribution by metal
    cur.execute("""
        SELECT metal, oxidation_state, COUNT(*) as cnt 
        FROM complexes 
        WHERE oxidation_state IS NOT NULL
        GROUP BY metal, oxidation_state 
        ORDER BY metal, oxidation_state
    """)
    ox_data = cur.fetchall()
    if ox_data:
        st.markdown("### Oxidation state distribution")
        for metal, ox, cnt in ox_data:
            st.markdown(f"**{metal}({ox}):** {cnt} complexes")


# ==================== STRUCTURES PAGE ====================
elif page == "⚗️ Structures":
    st.title("⚗️ Structures")
    st.markdown("Browse and download structural files for coordination compounds.")
    
    conn = get_db()
    cur = conn.cursor()
    
    # Filter
    metal_filter = st.selectbox("Filter by metal", ["All", "Ru", "Ir", "Rh", "Os", "Re"])
    
    if metal_filter == "All":
        cur.execute("SELECT id, metal, oxidation_state, smiles_ligands, metal_smiles, selfies FROM complexes ORDER BY id LIMIT 50")
    else:
        cur.execute("SELECT id, metal, oxidation_state, smiles_ligands, metal_smiles, selfies FROM complexes WHERE metal = ? ORDER BY id LIMIT 50", (metal_filter,))
    
    complexes = cur.fetchall()
    
    st.markdown(f"### Showing {len(complexes)} complexes")
    
    for cid, metal, ox, smiles, metal_smiles, selfies in complexes:
        with st.expander(f"#{cid} — {metal}({ox})"):
            col1, col2 = st.columns([1, 2])
            
            with col1:
                img = render_mol_image(metal_smiles or smiles, size=(300, 250))
                if img:
                    st.image(img)
            
            with col2:
                # Show donor atoms
                cur2 = conn.cursor()
                cur2.execute("SELECT donor_atoms FROM complexes WHERE id = ?", (cid,))
                da_result = cur2.fetchone()
                donor_atoms = da_result[0] if da_result else None
                
                st.markdown(f"**MetalCytoSMILES:**")
                st.code(metal_smiles or smiles, language="text")
                
                if donor_atoms and donor_atoms != "None":
                    st.markdown(f"**Donor atoms:** {donor_atoms}")
                
                if selfies:
                    st.markdown(f"**SELFIES:**")
                    st.code(selfies[:200], language="text")
                
                # Downloads
                mol_path = os.path.join(MOL_DIR, f"complex_{cid}.mol")
                if os.path.exists(mol_path):
                    with open(mol_path, 'r') as f:
                        st.download_button("📥 MOL (V2000)", f.read(), f"complex_{cid}.mol", "chemical/x-mdl-molfile", key=f"sv_{cid}")
                
                mol3_path = os.path.join(MOL3_DIR, f"complex_{cid}.mol")
                if os.path.exists(mol3_path):
                    with open(mol3_path, 'r') as f:
                        st.download_button("📥 MOL (V3000)", f.read(), f"complex_{cid}_v3000.mol", "chemical/x-mdl-molfile", key=f"s3_{cid}")
