import os
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

# Paths - Linux version
project_root = Path("/home/odalo/Potato_Pathogen")
fasta_file = project_root / "data/raw/16S/GSR-DB_V3-V4_cluster-1_seqs_dada2.fasta"
output_dir = project_root / "results/16S"
figures_dir = project_root / "results/figures"

# Create directories
output_dir.mkdir(parents=True, exist_ok=True)
figures_dir.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("16S DADA2 File Analysis")
print("=" * 60)

# Step 1: Verify file location
print("\n1. Verifying file location...")
if fasta_file.exists():
    size_mb = fasta_file.stat().st_size / (1024 * 1024)
    print(f"   ✓ File found: {fasta_file}")
    print(f"   Size: {size_mb:.2f} MB")
else:
    print(f"   ✗ File not found: {fasta_file}")
    exit(1)

# Step 2: Load sequences
print("\n2. Loading sequences...")
sequences = list(SeqIO.parse(fasta_file, "fasta"))
print(f"   Total sequences: {len(sequences)}")

# Step 3: Sequence statistics
print("\n3. Calculating sequence statistics...")
sequence_lengths = [len(seq.seq) for seq in sequences]
mean_length = sum(sequence_lengths) / len(sequence_lengths)
min_length = min(sequence_lengths)
max_length = max(sequence_lengths)

print(f"   Mean length: {mean_length:.0f} bp")
print(f"   Min length: {min_length} bp")
print(f"   Max length: {max_length} bp")

# Create sequence stats dataframe
sequence_data = []
for seq in sequences:
    sequence_data.append({
        'sequence_id': seq.id,
        'length': len(seq.seq),
        'description': seq.description[:100]
    })

df = pd.DataFrame(sequence_data)
df.to_csv(output_dir / "16s_sequence_stats.csv", index=False)
print(f"   ✓ Saved: {output_dir / '16s_sequence_stats.csv'}")

# Step 4: Extract amplicon region
print("\n4. Analyzing amplicon regions...")
def extract_sequence_type(description):
    if 'V3' in description and 'V4' in description:
        return 'V3-V4'
    elif 'V3' in description:
        return 'V3'
    elif 'V4' in description:
        return 'V4'
    else:
        return 'Unknown'

region_counts = Counter()
for seq in sequences:
    region = extract_sequence_type(seq.description)
    region_counts[region] += 1

print("   Amplicon region distribution:")
for region, count in region_counts.items():
    print(f"     {region}: {count} sequences")

# Step 5: Length distribution plot
print("\n5. Generating length distribution plot...")
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(sequence_lengths, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
ax.set_xlabel('Sequence Length (bp)')
ax.set_ylabel('Frequency')
ax.set_title('16S DADA2 Sequence Length Distribution')
ax.axvline(x=mean_length, color='red', linestyle='--', label=f'Mean: {mean_length:.0f} bp')
ax.legend()
plt.tight_layout()
plt.savefig(figures_dir / "16s_length_distribution.png", dpi=150)
plt.close()
print(f"   ✓ Saved: {figures_dir / '16s_length_distribution.png'}")

# Step 6: GC content analysis
print("\n6. Calculating GC content...")
def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

gc_contents = [calculate_gc_content(seq.seq) for seq in sequences]
mean_gc = sum(gc_contents) / len(gc_contents)

print(f"   Mean GC content: {mean_gc:.1f}%")
print(f"   Min GC content: {min(gc_contents):.1f}%")
print(f"   Max GC content: {max(gc_contents):.1f}%")

# GC content plot
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(gc_contents, bins=50, edgecolor='black', alpha=0.7, color='green')
ax.set_xlabel('GC Content (%)')
ax.set_ylabel('Frequency')
ax.set_title('16S Sequence GC Content Distribution')
ax.axvline(x=mean_gc, color='red', linestyle='--', label=f'Mean: {mean_gc:.1f}%')
ax.legend()
plt.tight_layout()
plt.savefig(figures_dir / "16s_gc_content.png", dpi=150)
plt.close()
print(f"   ✓ Saved: {figures_dir / '16s_gc_content.png'}")

# Step 7: Categorize by length
print("\n7. Categorizing sequences by length...")
short_seqs = [seq for seq in sequences if len(seq.seq) < 400]
medium_seqs = [seq for seq in sequences if 400 <= len(seq.seq) <= 450]
long_seqs = [seq for seq in sequences if len(seq.seq) > 450]

print(f"   Short sequences (<400bp): {len(short_seqs)}")
print(f"   Medium sequences (400-450bp): {len(medium_seqs)}")
print(f"   Long sequences (>450bp): {len(long_seqs)}")

# Save categorized sequences
if short_seqs:
    SeqIO.write(short_seqs, output_dir / "16s_short_sequences.fasta", "fasta")
if medium_seqs:
    SeqIO.write(medium_seqs, output_dir / "16s_medium_sequences.fasta", "fasta")
if long_seqs:
    SeqIO.write(long_seqs, output_dir / "16s_long_sequences.fasta", "fasta")

# Step 8: Create summary
print("\n8. Creating analysis summary...")
summary_data = {
    'Metric': [
        'Total sequences',
        'Mean length (bp)',
        'Min length (bp)',
        'Max length (bp)',
        'Mean GC content (%)',
        'V3-V4 region sequences',
        'Short sequences (<400bp)',
        'Medium sequences (400-450bp)',
        'Long sequences (>450bp)'
    ],
    'Value': [
        len(sequences),
        f"{mean_length:.0f}",
        min_length,
        max_length,
        f"{mean_gc:.1f}",
        region_counts.get('V3-V4', 0),
        len(short_seqs),
        len(medium_seqs),
        len(long_seqs)
    ]
}

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(output_dir / "16s_analysis_summary.csv", index=False)
print("\n   16S Analysis Summary:")
print("   " + "=" * 40)
for _, row in summary_df.iterrows():
    print(f"   {row['Metric']}: {row['Value']}")
print("   " + "=" * 40)

print("\n" + "=" * 60)
print("16S Analysis Complete!")
print("=" * 60)
