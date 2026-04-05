import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
#REMEMBER TO UPDATE THE PATHS TO YOUR OWN SYSTEM

# Paths
project_root = Path("/home/odalo/Potato_Pathogen")
feature_path = project_root / "data/processed/ITS/feature_table.tsv"
taxonomy_path = project_root / "data/processed/ITS/taxonomy_ITS.csv"
results_dir = project_root / "results"

(results_dir / "figures").mkdir(parents=True, exist_ok=True)

print("Taxonomy Bar Plot")
print("-" * 60)

# 1. Load data
print("\n1. Loading data...")
feature_table = pd.read_csv(feature_path, sep='\t', skiprows=1, index_col=0).T
taxonomy = pd.read_csv(taxonomy_path, sep='\t')

print(f"   Feature table: {feature_table.shape[0]} samples x {feature_table.shape[1]} features")
print(f"   Taxonomy: {taxonomy.shape[0]} features")
print(f"   Taxonomy columns: {taxonomy.columns.tolist()}")

taxonomy.columns = ['FeatureID', 'Taxon', 'Confidence']
print(f"   Renamed columns: {taxonomy.columns.tolist()}")

# 2. Extract phylum
print("\n2. Extracting phylum from taxonomy...")
taxonomy['Phylum'] = taxonomy['Taxon'].str.extract(r'p__([^;]+)')
print(f"   Found {taxonomy['Phylum'].nunique()} unique phyla")

taxonomy_clean = taxonomy.dropna(subset=['Phylum'])
print(f"   Features with phylum assignment: {len(taxonomy_clean)}")

common_features = pd.Index(taxonomy_clean['FeatureID']).intersection(feature_table.columns)
print(f"   Common features between tables: {len(common_features)}")

if len(common_features) == 0:
    print("\n   Checking feature IDs...")
    print(f"   First 5 feature table columns: {feature_table.columns[:5].tolist()}")
    print(f"   First 5 taxonomy FeatureIDs: {taxonomy_clean['FeatureID'].head().tolist()}")
else:
    feature_table_filtered = feature_table[common_features]
    taxonomy_filtered = taxonomy_clean.set_index('FeatureID').loc[common_features]

    # 3. Aggregate by phylum
    print("\n3. Aggregating by phylum...")
    feature_by_phylum = feature_table_filtered.T.groupby(taxonomy_filtered['Phylum']).sum().T
    relative_abundance = feature_by_phylum.div(feature_by_phylum.sum(axis=1), axis=0) * 100

    top_phyla = relative_abundance.sum().sort_values(ascending=False).head(10).index
    relative_abundance_top = relative_abundance[top_phyla].copy()

    other_cols = relative_abundance.columns[~relative_abundance.columns.isin(top_phyla)]
    relative_abundance_top['Other'] = relative_abundance[other_cols].sum(axis=1)

    print("\n   Top 10 phyla by abundance:")
    for phylum in top_phyla:
        print(f"     - {phylum}: {relative_abundance[phylum].mean():.1f}%")

    # 4. Bar plot
    print("\n4. Generating bar plot...")
    fig, ax = plt.subplots(figsize=(14, 8))
    relative_abundance_top.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
    ax.set_xlabel('Sample', fontsize=12)
    ax.set_ylabel('Relative Abundance (%)', fontsize=12)
    ax.set_title('Fungal Community Composition at Phylum Level', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(results_dir / "figures/taxonomy_barplot.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Saved: {results_dir / 'figures/taxonomy_barplot.png'}")

    # 5. Genus heatmap
    print("\n5. Creating genus-level heatmap...")
    taxonomy['Genus'] = taxonomy['Taxon'].str.extract(r'g__([^;]+)')
    taxonomy_genus = taxonomy.dropna(subset=['Genus'])

    common_genus = pd.Index(taxonomy_genus['FeatureID']).intersection(feature_table.columns)

    if len(common_genus) > 0:
        feature_table_genus = feature_table[common_genus]
        taxonomy_genus_filtered = taxonomy_genus.set_index('FeatureID').loc[common_genus]

        feature_by_genus = feature_table_genus.T.groupby(taxonomy_genus_filtered['Genus']).sum().T
        relative_genus = feature_by_genus.div(feature_by_genus.sum(axis=1), axis=0) * 100

        top_genera = relative_genus.sum().sort_values(ascending=False).head(20).index
        relative_genus_top = relative_genus[top_genera]

        fig, ax = plt.subplots(figsize=(12, 10))
        im = ax.imshow(relative_genus_top.T, aspect='auto', cmap='viridis')
        ax.set_xticks(range(len(relative_genus_top.index)))
        ax.set_xticklabels(relative_genus_top.index, rotation=45, ha='right', fontsize=8)
        ax.set_yticks(range(len(relative_genus_top.columns)))
        ax.set_yticklabels(relative_genus_top.columns, fontsize=9)
        ax.set_xlabel('Samples', fontsize=12)
        ax.set_ylabel('Genus', fontsize=12)
        ax.set_title('Top 20 Genera Relative Abundance (%)', fontsize=14)
        plt.colorbar(im, ax=ax, label='Relative Abundance (%)')
        plt.tight_layout()
        plt.savefig(results_dir / "figures/genus_heatmap.png", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"   Saved: {results_dir / 'figures/genus_heatmap.png'}")
    else:
        print("   No genus data available.")

print("\nTaxonomy analysis complete.")
print("-" * 60)