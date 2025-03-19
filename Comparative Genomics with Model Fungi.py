#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
import requests
import gzip
import shutil
import subprocess
import concurrent.futures
from tqdm import tqdm
import time

# Configuration
output_dir = "comparative_analysis"
model_species_dir = os.path.join(output_dir, "model_fungi")
local_blast_db = os.path.join(output_dir, "blast_db")
max_workers = 4  # Adjust based on CPU cores

# Model fungi with well-annotated genomes (including species that form mycorrhizal associations)
model_fungi = {
    "Rhizophagus irregularis": "GCF_000439145.1",  # Arbuscular mycorrhizal fungus
    "Laccaria bicolor": "GCF_000143565.1",  # Ectomycorrhizal fungus
    "Tuber melanosporum": "GCF_000151645.1",  # Truffle (ectomycorrhizal)
    "Glomus versiforme": "GCA_003550325.1",  # Arbuscular mycorrhizal fungus
    "Gigaspora rosea": "GCA_004153825.1",  # Arbuscular mycorrhizal fungus
    "Rhizoctonia solani": "GCF_000524565.1",  # Plant pathogen (for comparison)
    "Neurospora crassa": "GCF_000182925.2",  # Model organism
    "Saccharomyces cerevisiae": "GCF_000146045.2"  # Model yeast
}


def setup_environment():
    """Create output directories and check BLAST installation"""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(model_species_dir, exist_ok=True)
    os.makedirs(local_blast_db, exist_ok=True)

    # Check if BLAST is installed
    try:
        subprocess.run(["blastp", "-version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("BLAST is installed and available")
    except FileNotFoundError:
        print("BLAST not found in PATH. Please install BLAST+ tools")
        print("You can install it with: sudo apt-get install ncbi-blast+")
        return False

    return True


def download_reference_proteomes():
    """Download proteomes for model fungi"""
    print("Downloading reference proteomes for model fungi...")

    for species, assembly in model_fungi.items():
        species_dir = os.path.join(model_species_dir, assembly)
        os.makedirs(species_dir, exist_ok=True)

        # Download protein FASTA
        protein_file = os.path.join(species_dir, f"{assembly}_protein.faa.gz")
        protein_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{assembly[4:7]}/{assembly[7:10]}/{assembly[10:13]}/{assembly}/protein.faa.gz"

        if not os.path.exists(protein_file):
            try:
                print(f"Downloading {species} ({assembly}) proteome...")
                response = requests.get(protein_url, stream=True)
                if response.status_code == 404:
                    # Try GCA if GCF not found
                    alt_id = f"GCA_{assembly[4:]}"
                    protein_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{alt_id[4:7]}/{alt_id[7:10]}/{alt_id[10:13]}/{alt_id}/protein.faa.gz"
                    response = requests.get(protein_url, stream=True)

                with open(protein_file, 'wb') as f:
                    shutil.copyfileobj(response.raw, f)

                print(f"Downloaded {protein_file}")

                # Add a small delay to avoid overwhelming the server
                time.sleep(1)
            except Exception as e:
                print(f"Error downloading {species} proteome: {e}")

    # Extract and combine all model fungi proteins
    combined_fasta = os.path.join(model_species_dir, "all_model_fungi.faa")
    with open(combined_fasta, 'w') as out_f:
        for species, assembly in model_fungi.items():
            protein_file = os.path.join(model_species_dir, assembly, f"{assembly}_protein.faa.gz")

            if os.path.exists(protein_file):
                try:
                    with gzip.open(protein_file, 'rt') as in_f:
                        # Modify headers to include species name
                        for line in in_f:
                            if line.startswith('>'):
                                # Get the original protein ID
                                protein_id = line.strip().split()[0][1:]
                                # Create a new header with species
                                new_header = f">{species}|{protein_id} {' '.join(line.strip().split()[1:])}\n"
                                out_f.write(new_header)
                            else:
                                out_f.write(line)
                except Exception as e:
                    print(f"Error processing {protein_file}: {e}")

    print(f"Created combined reference proteome: {combined_fasta}")
    return combined_fasta


def create_blast_db(fasta_file):
    """Create BLAST database from FASTA file"""
    db_name = os.path.join(local_blast_db, os.path.basename(fasta_file).split('.')[0])

    cmd = f"makeblastdb -in {fasta_file} -dbtype prot -out {db_name}"

    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"Created BLAST database: {db_name}")
        return db_name
    except subprocess.CalledProcessError as e:
        print(f"Error creating BLAST database: {e}")
        return None


def run_blast(query_fasta, db_name, output_file, evalue=1e-5):
    """Run BLAST search against the model fungi database"""
    cmd = NcbiblastpCommandline(
        query=query_fasta,
        db=db_name,
        out=output_file,
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
        evalue=evalue,
        max_target_seqs=5
    )

    try:
        stdout, stderr = cmd()
        print(f"BLAST search completed: {output_file}")
        return output_file
    except Exception as e:
        print(f"Error running BLAST: {e}")
        return None


def find_receptor_candidates(assembly_dir, model_db, output_dir):
    """Find potential strigalactone receptor candidates using BLAST"""
    # Get query proteins
    protein_files = [f for f in os.listdir(assembly_dir) if f.endswith("_proteins.faa")]

    if not protein_files:
        print(f"No protein files found in {assembly_dir}")
        return

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Known strigalactone-related protein keywords
    strigalactone_keywords = [
        "strigolactone", "karrikin", "d14", "max2", "smax1", "smxl",
        "alpha/beta hydrolase", "f-box", "plant hormone", "carotenoid cleavage",
        "DWARF14", "KARRIKIN INSENSITIVE", "MORE AXILLARY GROWTH"
    ]

    all_candidates = []

    # Process each assembly
    for protein_file in protein_files:
        assembly_acc = protein_file.split("_proteins.faa")[0]
        query_file = os.path.join(assembly_dir, protein_file)

        # Skip if already processed
        result_file = os.path.join(output_dir, f"{assembly_acc}_blast_results.tsv")
        if os.path.exists(result_file):
            print(f"BLAST results for {assembly_acc} already exist, skipping")

            # Load existing results
            try:
                blast_df = pd.read_csv(result_file, sep='\t', header=None)
                blast_df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch",
                                    "gapopen", "qstart", "qend", "sstart", "send",
                                    "evalue", "bitscore", "stitle"]
            except Exception as e:
                print(f"Error loading existing results for {assembly_acc}: {e}")
                continue
        else:
            # Run BLAST
            blast_result = run_blast(query_file, model_db, result_file)
            if not blast_result:
                continue

            # Load BLAST results
            try:
                blast_df = pd.read_csv(result_file, sep='\t', header=None)
                blast_df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch",
                                    "gapopen", "qstart", "qend", "sstart", "send",
                                    "evalue", "bitscore", "stitle"]
            except Exception as e:
                print(f"Error loading BLAST results for {assembly_acc}: {e}")
                continue

        print(f"Analyzing BLAST results for {assembly_acc}")

        # Find potential strigalactone-related proteins based on BLAST hits
        candidates = []
        for _, row in blast_df.iterrows():
            hit_title = row["stitle"].lower() if isinstance(row["stitle"], str) else ""

            # Check if hit matches any strigalactone-related keywords
            if any(keyword.lower() in hit_title for keyword in strigalactone_keywords):
                candidates.append({
                    "assembly": assembly_acc,
                    "protein_id": row["qseqid"],
                    "hit_id": row["sseqid"],
                    "percent_identity": row["pident"],
                    "evalue": row["evalue"],
                    "bitscore": row["bitscore"],
                    "hit_description": row["stitle"],
                    "matching_keyword": next((kw for kw in strigalactone_keywords if kw.lower() in hit_title), "")
                })

        if candidates:
            # Add to all candidates
            all_candidates.extend(candidates)

            # Save candidates for this assembly
            candidates_df = pd.DataFrame(candidates)
            candidates_file = os.path.join(output_dir, f"{assembly_acc}_receptor_candidates.tsv")
            candidates_df.to_csv(candidates_file, sep='\t', index=False)
            print(f"Found {len(candidates)} potential receptor candidates for {assembly_acc}")
        else:
            print(f"No receptor candidates found for {assembly_acc}")

    # Save all candidates
    if all_candidates:
        all_df = pd.DataFrame(all_candidates)
        all_file = os.path.join(output_dir, "all_receptor_candidates.tsv")
        all_df.to_csv(all_file, sep='\t', index=False)
        print(f"Found a total of {len(all_df)} potential receptor candidates across all assemblies")
    else:
        print("No receptor candidates found")

    return all_candidates


def extract_candidate_sequences(all_candidates, input_dir, output_dir):
    """Extract sequences for receptor candidates"""
    if not all_candidates:
        print("No candidates to extract sequences for")
        return

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Create dictionary of protein_id to assembly
    protein_map = {candidate["protein_id"]: candidate["assembly"] for candidate in all_candidates}

    # Load all protein sequences
    all_proteins = {}
    for assembly in set(protein_map.values()):
        protein_file = os.path.join(input_dir, f"{assembly}_proteins.faa")

        if os.path.exists(protein_file):
            try:
                for record in SeqIO.parse(protein_file, "fasta"):
                    all_proteins[record.id] = record
            except Exception as e:
                print(f"Error loading proteins from {protein_file}: {e}")

    # Extract candidate sequences
    candidates_file = os.path.join(output_dir, "receptor_candidates.faa")
    with open(candidates_file, 'w') as f:
        for protein_id in protein_map.keys():
            if protein_id in all_proteins:
                record = all_proteins[protein_id]
                f.write(f">{record.id}\n{record.seq}\n")

    print(f"Extracted {len(protein_map)} candidate sequences to {candidates_file}")
    return candidates_file


if __name__ == "__main__":
    # Setup environment
    if not setup_environment():
        print("Environment setup failed")
        exit(1)

    # Download reference proteomes
    ref_proteome = download_reference_proteomes()

    # Create BLAST database
    db_name = create_blast_db(ref_proteome)

    if not db_name:
        print("Failed to create BLAST database")
        exit(1)

    # Find receptor candidates
    candidates = find_receptor_candidates("receptor_candidates", db_name, os.path.join(output_dir, "blast_results"))

    # Extract candidate sequences
    if candidates:
        extract_candidate_sequences(candidates, "receptor_candidates", os.path.join(output_dir, "structure_candidates"))