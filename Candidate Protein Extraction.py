#!/usr/bin/env python
import os
import pandas as pd
import gffutils
from Bio import SeqIO
import gzip
import concurrent.futures
from tqdm import tqdm

# Configuration
download_dir = "/Volumes/Elements/Assembly"
metadata_csv = "ncbi_amf_assembly_details.csv"
output_dir = "receptor_candidates"
max_workers = 4  # Adjust based on your CPU

# Create output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created directory: {output_dir}")

# We'll use a more general approach since the GFFs are poorly annotated
# Instead of relying heavily on keywords, we'll extract all proteins and compare them with model species

# Load metadata to identify assemblies with GFFs
metadata = pd.read_csv(metadata_csv)
valid_assemblies = []

for _, row in metadata.iterrows():
    assembly_acc = str(row["AssemblyAccession"]).strip()
    gff_filepath = os.path.join(download_dir, f"{assembly_acc}.gff.gz")
    fasta_filepath = os.path.join(download_dir, f"{assembly_acc}.fna.gz")

    # Check if both GFF and FASTA files exist
    if os.path.exists(gff_filepath) and os.path.exists(fasta_filepath):
        valid_assemblies.append({
            "assembly_acc": assembly_acc,
            "assembly_name": row["AssemblyName"],
            "gff_filepath": gff_filepath,
            "fasta_filepath": fasta_filepath
        })

print(f"Found {len(valid_assemblies)} assemblies with both GFF and FASTA files")


def extract_protein_sequences(assembly_info):
    """Extract all protein-coding sequences from the genome for comparative analysis"""
    assembly_acc = assembly_info["assembly_acc"]
    fasta_filepath = assembly_info["fasta_filepath"]
    gff_filepath = assembly_info["gff_filepath"]

    # Output files
    protein_fasta = os.path.join(output_dir, f"{assembly_acc}_proteins.faa")
    protein_metadata = os.path.join(output_dir, f"{assembly_acc}_protein_metadata.tsv")

    # Skip if already processed
    if os.path.exists(protein_fasta) and os.path.exists(protein_metadata):
        print(f"Assembly {assembly_acc} already processed, skipping")
        return

    # Load or create GFF database
    db_path = gff_filepath + ".db"
    try:
        if not os.path.exists(db_path):
            print(f"Creating GFF database for {assembly_acc}")
            db = gffutils.create_db(gff_filepath, dbfn=db_path, force=True,
                                    keep_order=True, merge_strategy="merge",
                                    sort_attribute_values=True)
        else:
            db = gffutils.FeatureDB(db_path, keep_order=True)
    except Exception as e:
        print(f"Error with GFF database for {assembly_acc}: {e}")
        return

    # Load FASTA genome
    try:
        with gzip.open(fasta_filepath, "rt") as handle:
            fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    except Exception as e:
        print(f"Error loading FASTA from {fasta_filepath}: {e}")
        return

    candidate_proteins = {}
    candidate_metadata = []

    # Extract all gene/mRNA features
    print(f"Extracting all protein-coding genes from {assembly_acc}")
    gene_count = 0

    # Process each gene in the database
    for gene in db.features_of_type('gene'):
        gene_count += 1
        gene_id = gene.id
        gene_description = ""

        # Try to get a description from gene attributes
        for key, values in gene.attributes.items():
            if key in ['product', 'description', 'Note', 'name']:
                for value in values:
                    if value and len(value) > 0:
                        gene_description = value
                        break

        # Get mRNAs associated with this gene
        mrnas = list(db.children(gene, featuretype='mRNA'))
        if not mrnas:
            # If no mRNAs, try to find CDS features directly under gene
            cds_features = list(db.children(gene, featuretype='CDS'))
            if cds_features:
                # Create a virtual mRNA
                mrna_id = f"{gene_id}_mRNA"
                mrnas = [(gene, mrna_id)]

        # Process each mRNA to get its CDS
        for mrna in mrnas:
            try:
                if isinstance(mrna, tuple):
                    # Handle virtual mRNA case
                    gene_feature, mrna_id = mrna
                    cds_parts = list(db.children(gene_feature, featuretype='CDS'))
                else:
                    mrna_id = mrna.id
                    cds_parts = list(db.children(mrna, featuretype='CDS'))

                if not cds_parts:
                    continue

                # Concatenate CDS parts and translate
                cds_sequences = []
                for cds in sorted(cds_parts, key=lambda x: x.start):
                    if cds.seqid in fasta_dict:
                        # Get the sequence, accounting for strand
                        seq = fasta_dict[cds.seqid].seq[cds.start - 1:cds.end]
                        if cds.strand == '-':
                            seq = seq.reverse_complement()
                        cds_sequences.append(str(seq))

                if not cds_sequences:
                    continue

                # Concatenate and add to candidates
                complete_cds = ''.join(cds_sequences)

                # Skip sequences that aren't multiples of 3 (likely incomplete)
                if len(complete_cds) % 3 != 0:
                    continue

                protein_id = f"{assembly_acc}|{mrna_id}"

                # Add metadata
                candidate_metadata.append({
                    "protein_id": protein_id,
                    "assembly": assembly_acc,
                    "gene_id": gene_id,
                    "mrna_id": mrna_id,
                    "description": gene_description,
                    "cds_length": len(complete_cds)
                })

                # Store DNA sequence for later translation
                candidate_proteins[protein_id] = complete_cds
            except Exception as e:
                print(f"Error processing mRNA {mrna_id if 'mrna_id' in locals() else 'unknown'}: {e}")
                continue

    # Write proteins to FASTA
    with open(protein_fasta, 'w') as f_out:
        for protein_id, cds_seq in candidate_proteins.items():
            f_out.write(f">{protein_id}\n{cds_seq}\n")

    # Write metadata
    if candidate_metadata:
        df = pd.DataFrame(candidate_metadata)
        df.to_csv(protein_metadata, sep='\t', index=False)

    print(f"Processed {assembly_acc}: Found {len(candidate_proteins)} receptor candidates")
    return


# Process assemblies in parallel
with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
    list(tqdm(executor.map(extract_protein_sequences, valid_assemblies),
              total=len(valid_assemblies),
              desc="Processing assemblies"))

# Combine results
all_candidates = []
for assembly_info in valid_assemblies:
    assembly_acc = assembly_info["assembly_acc"]
    metadata_file = os.path.join(output_dir, f"{assembly_acc}_receptor_metadata.tsv")
    if os.path.exists(metadata_file):
        try:
            df = pd.read_csv(metadata_file, sep='\t')
            all_candidates.append(df)
        except Exception as e:
            print(f"Error reading {metadata_file}: {e}")

if all_candidates:
    combined_df = pd.concat(all_candidates, ignore_index=True)
    combined_df.to_csv(os.path.join(output_dir, "all_receptor_candidates.tsv"), sep='\t', index=False)
    print(f"Found a total of {len(combined_df)} receptor candidates across all assemblies")
else:
    print("No receptor candidates found")