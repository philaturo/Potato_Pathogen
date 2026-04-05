import pandas as pd
import zipfile
import tempfile
from pathlib import Path
#REMEMBER TO UPDATE THE PATHS TO YOUR OWN SYSTEM

def extract_qza_table(qza_path, output_csv):
    print(f"Extracting: {qza_path}")
    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(qza_path, 'r') as zip_ref:
            zip_ref.extractall(tmpdir)

        data_file = Path(tmpdir) / "data" / "feature-table.biom"
        if not data_file.exists():
            data_file = Path(tmpdir) / "data" / "data.biom"

        if data_file.exists():
            import biom
            table = biom.load_table(str(data_file))
            df = table.to_dataframe().T
            df.to_csv(output_csv, sep='\t')
            print(f"  Saved: {output_csv} ({df.shape[0]} samples x {df.shape[1]} features)")
            return df
        else:
            print(f"  ERROR: Could not find data file in {tmpdir}")
            return None


def extract_qza_taxonomy(qza_path, output_csv):
    print(f"Extracting: {qza_path}")
    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(qza_path, 'r') as zip_ref:
            zip_ref.extractall(tmpdir)

        tax_file = Path(tmpdir) / "data" / "taxonomy.tsv"
        if tax_file.exists():
            df = pd.read_csv(tax_file, sep='\t')
            df.to_csv(output_csv, sep='\t', index=False)
            print(f"  Saved: {output_csv} ({df.shape[0]} features)")
            return df
        else:
            print(f"  ERROR: Could not find taxonomy file")
            return None


# Paths
project_root = Path("/home/odalo/Potato_Pathogen")
qza_dir = Path("/home/odalo/Documents/bioinfo/reference/UNITE_ITS/sh_qiime_release_04.04.2024/developer")
output_dir = project_root / "data/processed/ITS"
output_dir.mkdir(parents=True, exist_ok=True)

print("Extracting QZA Files")
print("-" * 50)

table_qza = qza_dir / "table.qza"
if table_qza.exists():
    extract_qza_table(table_qza, output_dir / "feature_table.csv")
else:
    print(f"ERROR: {table_qza} not found")

tax_qza = qza_dir / "taxonomy.qza"
if tax_qza.exists():
    extract_qza_taxonomy(tax_qza, output_dir / "taxonomy_ITS.csv")
else:
    print(f"ERROR: {tax_qza} not found")

print("\nExtraction complete.")
print(f"Files saved in: {output_dir}")