import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set seaborn style
sns.set(style="whitegrid")

# Adjusting the species names to match your desired formatting
def format_species_name(species):
    # Replace underscores with spaces and ensure correct formatting for genus and species
    species = species.replace("_", " ")
    species_parts = species.split(" ")
    if len(species_parts) == 2:
        species = f"{species_parts[0]} {species_parts[1]}"
    return species

# Your reference genes
genes_of_interest = [
    'Helicase', 'RHIP1', 'UPL7',
    'F-Box', 'AP2M', 'YLS8', 'PP2AA3','EXPRS', 'MON1', 
    'His3.3', 'PP2A-1', 'UBC21','AT5G12240','AT4G33380'
]

# Load all species data
species_data = {
    'Malus_domestica': pd.read_csv('Malus_domestica_reference2.csv', index_col=0),
    'Fragaria_vesca': pd.read_csv('Fragaria_vesca_references2.csv', index_col=0),
    'Ulmus_minor': pd.read_csv('Ulmus_minor_references2.csv', index_col=0),
    'Hippophae_rhamnoides': pd.read_csv('Hippophae_rhamnoides_references2.csv', index_col=0),
    'Quercus_robur': pd.read_csv('Quercus_robur_references2.csv', index_col=0),
    'Castanea_mollissima': pd.read_csv('Castanea_mollissima_references2.csv', index_col=0),
    'Juglans_regia': pd.read_csv('Juglans_regia_references2.csv', index_col=0),
    'Carya_illinoinensis': pd.read_csv('Carya_illinoinensis_references2.csv', index_col=0),
    'Gynostemma_pentaphyllum': pd.read_csv('Gynostemma_pentaphyllum_references2.csv', index_col=0),
    'Momordica_charantia': pd.read_csv('Momordica_charantia2_references2.csv', index_col=0),
    'Luffa_cylindrica': pd.read_csv('Luffa_cylindrica1_references2.csv', index_col=0),
    'Cucumis_sativus': pd.read_csv('Cucumis_sativus_references2.csv', index_col=0),
    'Citrullus_lanatus': pd.read_csv('Citrullus_lanatus_references2.csv', index_col=0),
    'Benincasa_hispida': pd.read_csv('Benincasa_hispida1_references2.csv', index_col=0),
    'Cucurbita_pepo': pd.read_csv('Cucurbita_pepo.tpms_references2.csv', index_col=0),
    'Cucurbita_moschata':pd.read_csv('Cucurbita_moschata1_references.csv', index_col=0),
}


# Summarize expression (adding isoforms together)
species_summed_tpms = {}
per_species_cvs = []

for species, df in species_data.items():
    df = df[df.index.isin(genes_of_interest)]  # Keep only reference genes

    if df.empty:
        continue

    # Sum isoforms (group by identical gene names)
    df_summed = df.groupby(df.index).sum()

    species_summed_tpms[species] = df_summed

    # CV per gene
    cvs = df_summed.std(axis=1) / df_summed.mean(axis=1)
    per_species_cvs.append(cvs)

# Merge CVs across species
all_cvs = pd.concat(per_species_cvs, axis=1)
mean_cvs = all_cvs.mean(axis=1)

# Top 10 most stable reference genes
top_10_genes = mean_cvs.nsmallest(10).index.tolist()
print("Top 10 most stable reference genes:", top_10_genes)

# Set colors
palette = sns.color_palette("tab10", n_colors=10)
gene_colors = dict(zip(top_10_genes, palette))

# ---- PLOTTING ----
sns.set(style="whitegrid")

# Create figure with GridSpec: 4×5 cells, last column reserved for legend
fig = plt.figure(figsize=(24, 20))
gs = fig.add_gridspec(4, 5, width_ratios=[1,1,1,1,0.3], wspace=0.3, hspace=0.4)

nrows, ncols = 8, 2
# Flattened list of the 16 plotting axes
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols, figsize=(16, 23), sharey='row'
)

plt.subplots_adjust(
    left=0.08,
    right=0.80,   # space for legend
    top=0.95,
    bottom=0.18   # space for figure label + caption
)
#axes = [fig.add_subplot(gs[i, j]) for i in range(4) for j in range(4)]

# Flatten axes for easy iteration
flat_axes = axes.flatten()

for ax, (species, df) in zip(flat_axes, species_summed_tpms.items()):
    # select only your top genes
    df = df.loc[df.index.intersection(top_10_genes)]
    if df.empty:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=13)
    else:
        df = df.T
        df_rel = df / df.mean()
        df_rel['Sample'] = df_rel.index
        df_long = df_rel.melt(
            id_vars='Sample', var_name='Gene', value_name='Relative TPM'
        )
        sns.lineplot(
            data=df_long, x='Sample', y='Relative TPM',
            hue='Gene', palette=gene_colors, ax=ax,
            legend=False, linewidth=1
        )
    n = df_rel.shape[0] if not df.empty else 0
    species_name = species.replace('_', r'\ ')
    ax.set_title(rf"$\it{{{species_name}}}$ (n={n})", fontsize=14)
    ax.set_xlabel("")
    ax.set_ylabel("Relative TPM", fontsize=12)
    ax.set_ylim(0, 5)
    row_idx = flat_axes.tolist().index(ax) // ncols
    if row_idx == nrows - 1:
        ax.set_xlabel("Samples", fontsize=12)
    else:
        ax.set_xlabel("")

    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.grid(axis='x', visible=False)
    ax.axhline(1, color='black', linestyle='--', linewidth=2)

# remove any extra axes (if fewer than 16 species)
for extra_ax in flat_axes[len(species_summed_tpms):]:
    fig.delaxes(extra_ax)

from matplotlib.lines import Line2D

# Legend handles
handles = [Line2D([], [], color=col, lw=3) for gene, col in gene_colors.items()]
labels = list(gene_colors.keys())

# Create separate axis for the legend (x0, y0, width, height)
legend_ax = fig.add_axes([0.82, 0.2, 0.16,0.55])  # Adjust position as needed
legend_ax.axis('off')
legend_ax.legend(
    handles, labels, loc='center left', fontsize=11, frameon=True, handlelength=2.5, borderpad=0.8,labelspacing=0.8)

plt.tight_layout()
plt.savefig("Expression_top10_reference_genes.svg", dpi=300, bbox_inches='tight')
plt.show()
