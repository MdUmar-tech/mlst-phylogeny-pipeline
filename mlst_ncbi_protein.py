import os
import zipfile
import shutil
from pathlib import Path
from Bio import SeqIO

# Define MLST genes
mlst_genes = {"atpD", "gyrB", "recA", "rpoB", "trpB"}

zip_folder = "Selected_genomes_based_on_16s_rRNA"
output_folder = "mlst_extracted_proteins"

os.makedirs(output_folder, exist_ok=True)

for zip_file in Path(zip_folder).glob("*.zip"):
    strain_name = zip_file.stem
    strain_output_folder = os.path.join(output_folder, strain_name)
    os.makedirs(strain_output_folder, exist_ok=True)

    print(f"\n🧬 Processing strain: {strain_name}")

    found_genes = {gene: False for gene in mlst_genes}

    temp_extract_path = os.path.join("temp_extract", strain_name)
    os.makedirs(temp_extract_path, exist_ok=True)

    try:
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            zip_ref.extractall(temp_extract_path)
    except zipfile.BadZipFile:
        print("❌ Invalid ZIP file")
        shutil.rmtree(temp_extract_path)
        continue

    gbff_files = (
        list(Path(temp_extract_path).glob("ncbi_dataset/data/GCA*/*genomic.gbff")) +
        list(Path(temp_extract_path).glob("ncbi_dataset/data/GCF*/*genomic.gbff"))
    )

    if not gbff_files:
        print("⚠️ No genomic.gbff found")
        shutil.rmtree(temp_extract_path)
        continue

    for gbff_file in gbff_files:
        print(f"  🔍 Scanning: {gbff_file.name}")

        with open(gbff_file, "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS" and "gene" in feature.qualifiers:
                        gene_name = feature.qualifiers["gene"][0]

                        if gene_name in mlst_genes and not found_genes[gene_name]:
                            if "translation" not in feature.qualifiers:
                                continue

                            protein_seq = feature.qualifiers["translation"][0]
                            out_file = os.path.join(
                                strain_output_folder, f"{gene_name}.faa"
                            )

                            with open(out_file, "a") as out_faa:
                                out_faa.write(
                                    f">{record.id}_{gene_name}\n{protein_seq}\n"
                                )

                            found_genes[gene_name] = True

    print("  📊 Extraction summary:")
    for gene in sorted(mlst_genes):
        status = "✅" if found_genes[gene] else "❌"
        print(f"     {gene:<6} {status}")

    shutil.rmtree(temp_extract_path)

print("\n🎉 MLST protein extraction completed!")
