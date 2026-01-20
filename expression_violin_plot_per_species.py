import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Define genes and species order
genes_of_interest = ['Z-ISO', 'ZDS', 'CrtISO','LCYE', 'LCYB', 'CYP97A', 'BCH', 'VDE', 'ZEP', 'NXS']
species_order = [
    'Malus domestica', 'Fragaria vesca', 'Ulmus minor', 'Hippophae rhamnoides',
    'Quercus robur', 'Castanea mollissima', 'Juglans regia', 'Carya illinoinensis','',
    'Gynostemma pentaphyllum','Momordica charantia', 'Luffa aegyptiaca','Cucumis sativus', 
    'Citrullus lanatus', 'Benincasa hispida', 'Cucurbita pepo','Cucurbita moschata'
]

# Define colors
species_colors = ['#1f78b4']*9 + ['#d95f02']*8

# Path to input folder
input_folder = Path("normalized_count_tables/")  

gene_mapping = {
    'Malus domestica': {
        #'PSY':['MD09G1281300','MD09G1146800','MD17G1133400','MD11G1010500','MD03G1007200'],
        #'PDS':['MD04G1023800'],
        'Z-ISO':['MD04G1028100'],
        'ZDS':['MD04G1220900','MD12G1237300'],
        'CrtISO':['MD14G1044400','MD08G1243700'],
        'LCYE':['MD02G1083500'],
        'LCYB':['MD00G1049000'],
        'CYP97A':['MD04G1239800','MD12G1257000'],
        #'CYP97C':['MD14G1097800'],
        'BCH':['MD01G1208300','MD07G1278900','MD06G1089200'],
        'VDE':['MD12G1255200'],
        'ZEP':['MD02G1172400','MD15G1284500'],
        'NXS':['MD00G1132400','MD09G1253000','MD06G1045700'],
    },
    'Fragaria vesca': {
        #'PSY':['XP_004289567.1','XP_004303895.1','XP_004296190.1'],
        #'PDS':['XP_004296964.1'],
        'Z-ISO':['XP_004298028.1'],
        'ZDS':['XP_004302025.1'],
        'CrtISO':['XP_004303382.1','XP_004301557.1'],
        'LCYE':['XP_004287582.1'],
        'LCYB':['XP_004303607.1'],
        'CYP97A':['XP_004302135.1'],
        #'CYP97C':['XP_004306170.1'],
        'BCH':['XP_004308054.1','XP_004300813.1'],
        'VDE':['XP_004302125.1'],
        'ZEP':['XP_004306384.1'],
        'NXS':['XP_004300312.1'],
    },
        'Ulmus minor': {
        #'PSY':['Umino31201.1','Umino03969.2','Umino18462.1'],
        #'PDS':['Umino51850.1'],
        'Z-ISO':['Umino51417.1'],
        'ZDS':['Umino29712.1','Umino29677.1','Umino29693.1'],
        'CrtISO':['Umino49671.1','Umino30776.1'],
        'LCYE':['Umino33628.1'],
        'LCYB':['Umino49479.1'],
        'CYP97A':['Umino30099.1'],
        #'CYP97C':['Umino50476.1'],
        'BCH':['Umino13563.1','Umino44146.1'],
        'VDE':['Umino30042.1'],
        'ZEP':['Umino38793.1'],
        'NXS':['Umino04414.1'],
    },
        'Hippophae rhamnoides': {
        #'PSY':['Hiprha1gene03092','Hiprha1gene06782'],
        #'PDS':['Hiprha1gene25899'],
        'Z-ISO':['Hiprha1gene09348'],
        'ZDS':['Hiprha1gene16000'],
        'CrtISO':['Hiprha1gene17015','Hiprha1gene21459','Hiprha1gene16344'],
        'LCYE':['Hiprha1gene10883','Hiprha1gene03401'],
        'LCYB':['Hiprha1gene27597'],
        'CYP97A':['Hiprha1gene26199'],
        #'CYP97C':['Hiprha1gene11236'],
        'BCH':['Hiprha1gene05693'],
        'VDE':['Hiprha1gene26091'],
        'ZEP':['Hiprha1gene15592'],
        'NXS':['Hiprha1gene00108','Hiprha1gene22818'],
    },
        'Quercus robur': {
        #'PSY':['QR_050393453.1','QR_050414418.1','QR_050414419.1'],
        #'PDS':['QR_050436637.1'],
        'Z-ISO':['QR_050399978.1','QR_050399979.1'],
        'ZDS':['QR_050419023.1'],
        'CrtISO':['QR_050388536.1','QR_050416925.1','QR_050416927.1'],
        'LCYE':['QR_050414224.1','QR_050414223.1'],
        'LCYB':['QR_050389511.1'],
        'CYP97A':['QR_050386579.1','QR_050386575.1','QR_050386572.1'],
        #'CYP97C':['QR_050387337.1','QR_050387336.1','QR_050387339.1','QR_050387338.1'],
        'BCH':['QR_050397113.1','QR_050410144.1','QR_050410142.1','QR_050410143.1'],
        'VDE':['QR_050386520.1','QR_050386525.1','QR_050386529.1'],
        'ZEP':['QR_050435776.1'],
        'NXS':['QR_050413833.1'],
    },
        'Castanea mollissima': {
        #'PSY':['KAF3946865.1','KAF3965429.1','KAF3965427.1'],
        #'PDS':['KAF3955907.1'],
        'Z-ISO':['KAF3959576.1'],
        'ZDS':['KAF3965006.1'],
        'CrtISO':['KAF3952049.1','KAF3956821.1'],
        'LCYE':['KAF3962181.1'],
        'LCYB':['KAF3946185.1','KAF3951486.1'],
        'CYP97A':['KAF3974451.1'],
        #'CYP97C':['KAF3974031.1'],
        'BCH':['KAF3953245.1','KAF3952196.1'],
        'VDE':['KAF3975440.1'],
        'ZEP':['KAF3963484.1'],
        'NXS':['KAF3949419.1','KAF3949418.1'],
    },
        'Juglans regia': {
        #'PSY':['XP_018816875.1','XP_035548768.1','XP_035551554.1','XP_018834528.1','XP_018834529.1','XP_018807983.1'],
        #'PDS':['XP_018828238.2'],
        'Z-ISO':['XP_018827677.1'],
        'ZDS':['XP_018845846.1'],
        'CrtISO':['XP_018848812.2','XP_018811007.1','XP_018811008.1'],
        'LCYE':['XP_018846198.1'],
        'LCYB':['XP_018817188.1'],
        'CYP97A':['XP_018817739.1'],
        #'CYP97C':['XP_018813031.1','XP_035544445.1'],
        'BCH':['XP_018809669.1','XP_018838613.1','XP_018852545.2'],
        'VDE':['XP_018826961.1'],
        'ZEP':['XP_018844974.1'],
        'NXS':['XP_018815488.2'],
    },
        'Carya illinoinensis': {
        #'PSY':['XP_042970360.1','XP_042972151.1','XP_042972152.1','XP_042954461.1'],
        #'PDS':['XP_042987139.1'],
        'Z-ISO':['XP_042971009.1'],
        'ZDS':['XP_042989449.1'],
        'CrtISO':['XP_042940746.1','XP_042983847.1'],
        'LCYE':['XP_042976990.1'],
        'LCYB':['XP_042942035.1'],
        'CYP97A':['XP_042988577.1'],
        #'CYP97C':['XP_042940175.1'],
        'BCH':['XP_042982516.1','XP_042985248.1','XP_042947713.1'],
        'VDE':['XP_042987638.1','XP_042987639.1'],
        'ZEP':['XP_042974075.1'],
        'NXS':['XP_042977548.1'],
    },
        'Gynostemma pentaphyllum': {
        #'PSY':['GP00731.1','GP00731.2','GP07893.1','GP35711.2','GP35711.1'],
        #'PDS':['GP18868.1'],
        'Z-ISO':['GP18982.1'],
        'ZDS':['GP23721.1','GP23785.1'],
        'CrtISO':['GP08653.1','GP36474.1'],
        'LCYE':['GP11616.1'],
        'LCYB':['GP14085.1'],
        'CYP97A':['GP36538.1'],
        #'CYP97C':['GP14775.1'],
        'BCH':['GP35469.1','GP20560.1'],
        'VDE':['GP11374.1'],
        'ZEP':['GP39143.1','GP39222.1'],
        'NXS':['GP15472.1','GP14998.1'],
    },
        'Momordica charantia': {
        #'PSY':['Moc04g30190.1','Moc10g09450.1','Moc02g00590.1'],
        #'PDS':['Moc05g10180.1'],
        'Z-ISO':['Moc05g08940.1'],
        'ZDS':['Moc06g26840.1'],
        'CrtISO':['Moc05g29080.1','Moc10g05510.1'],
        'LCYE':['Moc08g45730.1'],
        'LCYB':['Moc01g17190.1'],
        'CYP97A':['Moc10g06010.1'],
        #'CYP97C':['Moc01g22790.1'],
        'BCH':['Moc10g01250.1','Moc06g38230.1'],
        'VDE':['Moc08g38440.1'],
        'ZEP':['Moc02g13110.1'],
        'NXS':['Moc08g12330.1'],
    },
        'Luffa aegyptiaca': {
        #'PSY':['Lcy09g012870.1','Lcy11g006810.1','Lcy07g019180.1'],
        #'PDS':['Lcy10g007990.1'],
        'Z-ISO':['Lcy10g009130.1'],
        'ZDS':['Lcy03g004530.1'],
        'CrtISO':['Lcy06g018940.1','Lcy11g011320.1'],
        'LCYE':['Lcy08g018900.1'],
        'LCYB':['Lcy04g005840.1'],
        'CYP97A':['Lcy11g010510.1'],
        #'CYP97C':['Lcy04g011510.1'],
        'BCH':['Lcy11g001300.1','Lcy02g013470.1'],
        'VDE':['Lcy08g010620.1'],
        'ZEP':['Lcy07g004710.1'],
        'NXS':['Lcy02g001400.1'],
    },
        'Cucumis sativus': {
        #'PSY':['NC_026659.2_cds_XP_031741913.1_19472','NC_026659.2_cds_XP_004148907.1_19473','NC_026658.2_cds_XP_004147748.1_15202','NC_026658.2_cds_XP_031740055.1_15900','NC_026658.2_cds_XP_004142164.1_15901'],
        #'PDS':['NC_026658.2_cds_XP_031740146.1_14607','NC_026658.2_cds_XP_031740147.1_14608'],
        'Z-ISO':['NC_026658.2_cds_XP_004152101.1_14686'],
        'ZDS':['NC_026655.2_cds_XP_004142522.1_4378'],
        'CrtISO':['NC_026657.2_cds_XP_004136535.1_12638','NC_026657.2_cds_XP_011650571.1_8953'],
        'LCYE':['NC_026656.2_cds_XP_004141172.1_4581'],
        'LCYB':['NC_026658.2_cds_XP_004150761.1_15859'],
        'CYP97A':['NC_026657.2_cds_XP_004133753.1_9047'],
        #'CYP97C':['NC_026660.2_cds_XP_004143287.1_25297'],
        'BCH':['NC_026657.2_cds_XP_004140758.1_10464','NC_026659.2_cds_XP_004143973.1_18096'],
        'VDE':['NC_026656.2_cds_NP_001292657.1_5607'],
        'ZEP':['NC_026656.2_cds_NP_001292713.1_6342'],
        'NXS':['NC_026657.2_cds_XP_031738357.1_10909','NC_026657.2_cds_XP_004145860.1_10910','NC_026657.2_cds_XP_031738358.1_10911'],
    },
        'Citrullus lanatus': {
        #'PSY':['Cla97C01G008760.1','Cla97C07G137500.2','Cla97C02G050140.1'],
        #'PDS':['Cla97C07G142100.2'],
        'Z-ISO':['Cla97C07G142740.2'],
        'ZDS':['Cla97C06G118930.2'],
        'CrtISO':['Cla97C10G200950.2','Cla97C10G190930.2'],
        'LCYE':['Cla97C11G208040.2'],
        'LCYB':['Cla97C04G070940.1'],
        'CYP97A':['Cla97C10G191610.1'],
        #'CYP97C':['Cla97C04G076340.2'],
        'BCH':['Cla97C05G090480.1','Cla97C01G002480.2'],
        'VDE':['Cla97C11G216330.1'],
        'ZEP':['Cla97C02G038200.1'],
        'NXS':['Cla97C05G093720.2'],
    },
        'Benincasa hispida': {
        #'PSY':['BH_039019638.1','BH_039019639.1','BH_039018430.1','BH_039045277.1'],
        #'PDS':['BH_039041409.1'],
        'Z-ISO':['BH_039026618.1'],
        'ZDS':['BH_039023482.1'],
        'CrtISO':['BH_039047825.1','BH_039047826.1','BH_039049488.1'],
        'LCYE':['BH_039035238.1'],
        'LCYB':['BH_039026592.1'],
        'CYP97A':['BH_039050399.1'],
        #'CYP97C':['BH_039026659.1'],
        'BCH':['BH_039021520.1','BH_039051224.1'],
        'VDE':['BH_039034565.1'],
        'ZEP':['BH_039045392.1'],
        'NXS':['BH_039031262.1'],
    },
        'Cucurbita pepo': {
        #'PSY':['XP_023545264.1','XP_023550914.1','XP_023550915.1','XP_023553556.1','XP_023513201.1'],
        #'PDS':['XP_023539981.1','XP_023527611.1','XP_023527612.1','XP_023527613.1'],
        'Z-ISO':['XP_023540701.1'],
        'ZDS':['XP_023529959.1'],
        'CrtISO':['XP_023539389.1','XP_023553449.1'],
        'LCYE':['XP_023538575.1','XP_023546899.1'],
        'LCYB':['XP_023517264.1','XP_023545468.1'],
        'CYP97A':['XP_023552789.1'],
        #'CYP97C':['XP_023518813.1'],
        'BCH':['XP_023541483.1','XP_023551817.1','XP_023550966.1','XP_023550967.1'],
        'VDE':['XP_023536777.1'],
        'ZEP':['XP_023533320.1'],
        'NXS':['XP_023527642.1','XP_023527643.1'],
    },
        'Cucurbita moschata': {
        #'PSY':['NW_019268518.1_cds_XP_022963683.1_2548','NW_019268523.1_cds_XP_022938558.1_10216','NW_019268523.1_cds_XP_022938559.1_10217','NW_019268525.1_cds_XP_022941766.1_12887','NW_019268580.1_cds_XP_022931536.1_42172','NW_019268526.1_cds_XP_022943564.1_14809'],
        #'PDS':['NW_019268536.1_cds_XP_022955514.1_24853','NW_019268536.1_cds_XP_022955516.1_24854','NW_019268602.1_cds_XP_022932732.1_43183','NW_019268536.1_cds_XP_022955517.1_24855'],
        'Z-ISO':['NW_019268602.1_cds_XP_022932715.1_43131'],
        'ZDS':['NW_019268519.1_cds_XP_022933489.1_4994'],
        'CrtISO':['NW_019268517.1_cds_XP_022939505.1_1525','NW_019268567.1_cds_XP_022929364.1_40197'],
        'LCYE':['NW_019268522.1_cds_XP_022937731.1_8295','NW_019268528.1_cds_XP_022946860.1_15978'],
        'LCYB':['NW_019268554.1_cds_XP_022925158.1_36622','NW_019268614.1_cds_XP_022932998.1_43446'],
        'CYP97A':['NW_019268567.1_cds_XP_022929344.1_40259'],
        #'CYP97C':['NW_019268544.1_cds_XP_022962840.1_31044'],
        'BCH':['NW_019268527.1_cds_XP_022944865.1_15802','NW_019268550.1_cds_XP_022922621.1_34188','NW_019268566.1_cds_XP_022929023.1_40056'],
        'VDE':['NW_019268522.1_cds_XP_022937596.1_8786'],
        'ZEP':['NW_019268529.1_cds_XP_022947548.1_17205','NW_019268529.1_cds_XP_022947549.1_17206'],
        'NXS':['NW_019268520.1_cds_XP_022935036.1_5603','NW_019268520.1_cds_XP_022935038.1_5604','NW_019268520.1_cds_XP_022935039.1_5605'],
    },
}

plot_data = []

for species in species_order:
    species_file = species.replace(' ', '_') + '.txt'
    species_path = input_folder / species_file

    if not species_path.exists():
        print(f"Warning: File not found for {species}")
        continue

    df = pd.read_csv(species_path, sep="\t", index_col=0)

    if species not in gene_mapping:
        print(f" Warning: No gene mapping available for {species}")
        continue

    mapping = gene_mapping[species]
    df_mapped = {}

    for enzyme, gene_ids in mapping.items():
        found_genes = [g for g in gene_ids if g in df.index]
        if not found_genes:
            continue
        df_mapped[enzyme] = df.loc[found_genes].sum(axis=0)  # Sum of homologs per sample

    if not df_mapped:
        print(f"Warning: No genes mapped for {species}")
        continue

    df_mapped = pd.DataFrame(df_mapped).T

    for gene in genes_of_interest:
        if gene not in df_mapped.index:
            continue

        for sample_name, expr_value in df_mapped.loc[gene].items():
            plot_data.append({
                'Species': species,
                'Gene': gene,
                'Expression': expr_value
            })

# sanity check
if not plot_data:
    raise RuntimeError("No data was loaded—check your input_folder & gene_mapping!")

# build DataFrame
plot_df = pd.DataFrame(plot_data)

# enforce ordering
plot_df['Species'] = pd.Categorical(
    plot_df['Species'], categories=species_order, ordered=True
)
plot_df['Gene'] = pd.Categorical(
    plot_df['Gene'], categories=genes_of_interest, ordered=True
)

# ─── PLOTTING ──────────────────────────────────────────────────────────────────
sns.set(style="white")
n_genes = len(genes_of_interest)

#plot_df['Expression'] = np.log10(plot_df['Expression'] + 1)

fig, axes = plt.subplots(
    n_genes, 1,
    figsize=(3.5, 1.2 * n_genes),
    sharex=True,
    gridspec_kw={'hspace': 0.04}  
)
plt.subplots_adjust(right=0.8)
if n_genes == 1:
    axes = [axes]

for ax, gene in zip(axes, genes_of_interest):
    sub = plot_df[plot_df['Gene'] == gene]
    sns.violinplot(
        data=sub,
        x='Species',
        y='Expression',
        hue='Species',
        palette=species_colors,
        dodge=False,
        inner="box",
        linewidth=0,
        ax=ax       
    )
        # Add inner boxplot manually
    sns.boxplot(
        data=sub,
        x='Species',
        y='Expression',
        palette=species_colors,
        dodge=False,
        width=0.1,
        linewidth=0.5,
        fliersize=0,  # hide outliers if desired
        ax=ax
    )
    ax.set_title("")
    ax.set_ylabel(rf"$\it{{{gene}}}$", fontsize=8, fontweight='bold', rotation=0, ha='left', va='center')
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.tick_params(axis='y', which='both', length=0,labelsize=4)  # Hide the tick marks
    ax.tick_params(axis='x', which='both', length=0)
    ax.yaxis.set_label_coords(1.09, 0.5)
    ax.set_xlabel("")
    upper_limit = sub['Expression'].quantile(0.98)
    ax.set_ylim(0, upper_limit * 1.5)
    ax.grid(axis='x', visible=False)

    for side in ['top', 'right', 'bottom', 'left']:
        ax.spines[side].set_color('black')
        ax.spines[side].set_linewidth(0.5)

for i, (tick, species) in enumerate(zip(axes[-1].get_xticks(), species_order)):
    if i == 8:  # skip the 9th position (index 8)
        continue
    n = plot_df[plot_df['Species'] == species]['Expression'].groupby(plot_df['Gene']).count().max()
    if pd.notna(n) and n > 0:
        axes[0].text(
            tick,
            axes[0].get_ylim()[1] * 1.05,
            f"n={int(n)}",
            ha='left',
            va='bottom',
            fontsize=5,
            rotation=45
        )

# Restore black spines for both the main and twin axes
for ax in [axes[0]]:
    for side in ['top', 'right', 'bottom', 'left']:
        ax.spines[side].set_color('black')
        ax.spines[side].set_linewidth(0.5)

# remove the small per‐subplot legends
for ax in axes:
    if ax.get_legend():
        ax.get_legend().remove()

    # Remove legend from each plot
    if ax.get_legend():
        ax.get_legend().remove()

# Rotate x-axis labels
axes[-1].set_xticklabels(
    [
        r"$\it{{" + label.get_text().replace(' ', r'\ ') + "}}$"
        for label in axes[-1].get_xticklabels()
    ],
    rotation=45,
    ha='right',
    #va='bottom',
    fontsize=6
)

plt.tight_layout()
plt.savefig("carotenoid_violin_plots_per_lineage10.svg", dpi=300)
plt.show()
