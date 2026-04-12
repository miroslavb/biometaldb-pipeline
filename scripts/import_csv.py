#!/usr/bin/env python3
"""Import MetalCytoToxDB CSV into SQLite database."""

import csv
import io
import sqlite3
import os
import sys

import requests

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
CSV_PATH = os.path.join(DATA_DIR, "MetalCytoToxDB.csv")
DB_PATH = os.path.join(DATA_DIR, "biometaldb.sqlite")
ZENODO_URL = "https://zenodo.org/api/records/17106822/files/MetalCytoToxDB.csv/content"


def download_csv():
    """Download CSV from Zenodo if not present."""
    if os.path.exists(CSV_PATH):
        print(f"CSV already exists: {CSV_PATH}")
        return
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Downloading from Zenodo...")
    r = requests.get(ZENODO_URL, timeout=120)
    r.raise_for_status()
    with open(CSV_PATH, "w", encoding="utf-8") as f:
        f.write(r.text)
    print(f"Downloaded: {len(r.text)/(1024*1024):.1f} MB")


def create_database():
    """Create SQLite database from CSV."""
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
    
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    
    # Read CSV
    with open(CSV_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    print(f"Loaded {len(rows)} rows from CSV")
    
    # Create tables
    cur.executescript("""
        CREATE TABLE IF NOT EXISTS complexes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            smiles_ligands TEXT NOT NULL,
            metal TEXT NOT NULL,
            oxidation_state INTEGER,
            charge_complex INTEGER,
            metal_smiles TEXT,
            selfies TEXT,
            donor_atoms TEXT,
            UNIQUE(smiles_ligands, metal, oxidation_state)
        );
        
        CREATE TABLE IF NOT EXISTS measurements (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            complex_id INTEGER NOT NULL,
            abbreviation TEXT,
            ic50_dark REAL,
            ic50_dark_se REAL,
            ic50_dark_value REAL,
            ic50_light REAL,
            ic50_light_se REAL,
            ic50_light_value REAL,
            excitation_wavelength REAL,
            irradiation_time REAL,
            irradiation_power REAL,
            cell_line TEXT,
            exposure_time_h REAL,
            doi TEXT,
            year INTEGER,
            ic50_cisplatin REAL,
            ic50_cisplatin_se REAL,
            ic50_cisplatin_value REAL,
            counterion TEXT,
            atomic_number INTEGER,
            valence_e INTEGER,
            FOREIGN KEY (complex_id) REFERENCES complexes(id)
        );
        
        CREATE TABLE IF NOT EXISTS papers (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            doi TEXT UNIQUE,
            year INTEGER,
            title TEXT,
            authors TEXT
        );
        
        CREATE INDEX IF NOT EXISTS idx_complexes_metal ON complexes(metal);
        CREATE INDEX IF NOT EXISTS idx_complexes_smiles ON complexes(smiles_ligands);
        CREATE INDEX IF NOT EXISTS idx_measurements_cell_line ON measurements(cell_line);
        CREATE INDEX IF NOT EXISTS idx_measurements_doi ON measurements(doi);
        CREATE INDEX IF NOT EXISTS idx_measurements_complex_id ON measurements(complex_id);
        CREATE INDEX IF NOT EXISTS idx_papers_doi ON papers(doi);
    """)
    
    # Insert complexes
    complex_map = {}  # (smiles, metal, ox) -> id
    for row in rows:
        key = (row['SMILES_Ligands'], row['Metal'], 
               int(float(row['Oxidation_state'])) if row['Oxidation_state'] else None)
        if key not in complex_map:
            cur.execute("""
                INSERT OR IGNORE INTO complexes (smiles_ligands, metal, oxidation_state, charge_complex)
                VALUES (?, ?, ?, ?)
            """, (row['SMILES_Ligands'], row['Metal'], key[2],
                  int(float(row['Charge_complex'])) if row['Charge_complex'] else None))
            complex_map[key] = cur.lastrowid
    
    print(f"Inserted {len(complex_map)} unique complexes")
    
    # Insert measurements
    for row in rows:
        key = (row['SMILES_Ligands'], row['Metal'],
               int(float(row['Oxidation_state'])) if row['Oxidation_state'] else None)
        complex_id = complex_map[key]
        
        def safe_float(v):
            try:
                return float(v) if v and v.strip() else None
            except ValueError:
                return None
        
        def safe_int(v):
            try:
                return int(float(v)) if v and v.strip() else None
            except ValueError:
                return None
        
        cur.execute("""
            INSERT INTO measurements (
                complex_id, abbreviation, ic50_dark, ic50_dark_se, ic50_dark_value,
                ic50_light, ic50_light_se, ic50_light_value,
                excitation_wavelength, irradiation_time, irradiation_power,
                cell_line, exposure_time_h, doi, year,
                ic50_cisplatin, ic50_cisplatin_se, ic50_cisplatin_value,
                counterion, atomic_number, valence_e
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            complex_id,
            row.get('Abbreviation_in_the_article', ''),
            safe_float(row.get('IC50_Dark(M*10^-6)')),
            safe_float(row.get('IC50_Dark_standard_error(M*10^-6)')),
            safe_float(row.get('IC50_Dark_value')),
            safe_float(row.get('IC50_Light(M*10^-6)')),
            safe_float(row.get('IC50_Light_standard_error(M*10^-6)')),
            safe_float(row.get('IC50_Light_value')),
            safe_float(row.get('Excitation_Wavelength(nm)')),
            safe_float(row.get('Irradiation_Time(minutes)')),
            safe_float(row.get('Irradiation_Power(W*m^-2)')),
            row.get('Cell_line', ''),
            safe_float(row.get('Time(h)')),
            row.get('DOI', ''),
            safe_int(row.get('Year')),
            safe_float(row.get('IC50_Cisplatin(M*10^-6)')),
            safe_float(row.get('IC50_Cisplatin_standard_error(M*10^-6)')),
            safe_float(row.get('IC50_Cisplatin_value')),
            row.get('Counterion', ''),
            safe_int(row.get('Atomic_number')),
            safe_int(row.get('Valence_e')),
        ))
    
    print(f"Inserted {len(rows)} measurements")
    
    # Insert papers
    dois = set()
    for row in rows:
        doi = row.get('DOI', '')
        if doi and doi not in dois:
            dois.add(doi)
            cur.execute("INSERT OR IGNORE INTO papers (doi, year) VALUES (?, ?)",
                       (doi, safe_int(row.get('Year'))))
    
    print(f"Inserted {len(dois)} papers")
    
    conn.commit()
    
    # Verify
    cur.execute("SELECT COUNT(*) FROM complexes")
    n_complexes = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM measurements")
    n_measurements = cur.fetchone()[0]
    cur.execute("SELECT metal, COUNT(*) FROM complexes GROUP BY metal ORDER BY COUNT(*) DESC")
    metal_dist = cur.fetchall()
    
    print(f"\nDatabase verification:")
    print(f"  Complexes: {n_complexes}")
    print(f"  Measurements: {n_measurements}")
    print(f"  Metal distribution:")
    for metal, count in metal_dist:
        print(f"    {metal}: {count}")
    
    conn.close()
    print(f"\nDatabase saved: {DB_PATH}")
    print(f"Database size: {os.path.getsize(DB_PATH)/(1024*1024):.1f} MB")


if __name__ == "__main__":
    download_csv()
    create_database()
    print("\n✓ Stage 1 complete")
