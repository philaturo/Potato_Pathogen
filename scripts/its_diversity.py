import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu
from sklearn.manifold import MDS
from pathlib import Path

#REMEMBER TO UPDATE THE PATHS TO YOUR OWN SYSTEM
# Paths
project_root = Path("/home/odalo/Potato_Pathogen")
feature_path = project_root / "data/processed/ITS/feature_table.tsv"
taxonomy_path = project_root / "data/processed/ITS/taxonomy_ITS.csv"
metadata_path = project_root / "data/metadata/metadata.csv"
results_dir = project_root / "results"

(results_dir / "diversity").mkdir(parents=True, exist_ok=True)
(results_dir / "figures").mkdir(parents=True, exist_ok=True)

print("ITS Diversity Analysis")
print("-" * 60)

# 1. Load data
print("\n1. Loading data...")
feature_table = pd.read_csv(feature_path, sep='\t', skiprows=1, index_col=0).T
print(f"   Feature table: {feature_table.shape[0]} samples x {feature_table.shape[1]} features")

taxonomy = pd.read_csv(taxonomy_path)
print(f"   Taxonomy: {len(taxonomy)} features")

metadata = pd.read_csv(metadata_path)
print(f"   Metadata: {len(metadata)} samples")

# 2. Alpha diversity
print("\n2. Calculating alpha diversity...")
rel_abund = feature_table.div(feature_table.sum(axis=1), axis=0)
shannon = -np.sum(rel_abund * np.log(rel_abund + 1e-10), axis=1)
simpson = 1 - np.sum(rel_abund ** 2, axis=1)
observed = (feature_table > 0).sum(axis=1)

alpha_df = pd.DataFrame({
    'SampleID': feature_table.index,
    'Shannon': shannon.values,
    'Simpson': simpson.values,
    'Observed_ASVs': observed.values
})
alpha_df = alpha_df.merge(metadata, left_on='SampleID', right_on='FileID', how='left')
alpha_df.to_csv(results_dir / "diversity/alpha_diversity.csv", index=False)
print("   Alpha diversity saved.")

fig, ax = plt.subplots(figsize=(8, 6))
sns.boxplot(data=alpha_df, x='Type', y='Shannon', hue='Type',
            palette=['steelblue', 'orange'], legend=False, ax=ax)
ax.set_title('Shannon Diversity: Symptomatic vs Asymptomatic')
ax.set_ylabel('Shannon Diversity Index')
plt.tight_layout()
plt.savefig(results_dir / "figures/alpha_diversity.png", dpi=150)
plt.close()
print("   Alpha diversity plot saved.")

# 3. Beta diversity
print("\n3. Calculating beta diversity...")
distance_matrix = pdist(feature_table.values, metric='braycurtis')
square_distance = squareform(distance_matrix)

mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
pcoa_coords = mds.fit_transform(square_distance)

pcoa_df = pd.DataFrame(pcoa_coords, columns=['PCo1', 'PCo2'])
pcoa_df['SampleID'] = feature_table.index
pcoa_df = pcoa_df.merge(metadata, left_on='SampleID', right_on='FileID', how='left')

fig, ax = plt.subplots(figsize=(10, 8))
for disease_type in pcoa_df['Type'].unique():
    subset = pcoa_df[pcoa_df['Type'] == disease_type]
    ax.scatter(subset['PCo1'], subset['PCo2'], label=disease_type,
               s=80, alpha=0.7, edgecolors='black', linewidth=0.5)
ax.legend(title='Disease Status', fontsize=10)
ax.set_xlabel('PCo1', fontsize=12)
ax.set_ylabel('PCo2', fontsize=12)
ax.set_title('Bray-Curtis PCoA - Fungal Communities', fontsize=14)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(results_dir / "figures/beta_diversity_pcoa.png", dpi=150)
plt.close()
print("   Beta diversity plot saved.")

# 4. Statistical test
symptomatic = alpha_df[alpha_df['Type'] == 'Symptomatic']['Shannon']
asymptomatic = alpha_df[alpha_df['Type'] == 'Asymptomatic']['Shannon']
stat, p_value = mannwhitneyu(symptomatic, asymptomatic)

print("\n4. Statistical comparison:")
print(f"   Symptomatic vs Asymptomatic Shannon diversity")
print(f"   Mann-Whitney U test p-value: {p_value:.4f}")

print("\nDiversity analysis complete.")
print("-" * 60)
print("\nSummary:")
print(f"   Total samples: {len(feature_table)}")
print(f"   Total features: {feature_table.shape[1]}")
print(f"   Mean Shannon diversity: {shannon.mean():.2f} +/- {shannon.std():.2f}")
print(f"   Shannon diversity range: {shannon.min():.2f} - {shannon.max():.2f}")
print(f"\n   Symptomatic (n={len(symptomatic)}): {symptomatic.mean():.2f} +/- {symptomatic.std():.2f}")
print(f"   Asymptomatic (n={len(asymptomatic)}): {asymptomatic.mean():.2f} +/- {asymptomatic.std():.2f}")