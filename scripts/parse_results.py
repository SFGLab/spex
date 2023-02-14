import pandas as pd
from statannot import add_stat_annotation

# import matplotlib
import matplotlib.pyplot as plt
# import seaborn
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter

df = pd.read_csv("pcHi-C results SpEx - Log.csv", index_col=None)[["Note", "spearman"]]

#df = df.pivot_table(values='spearman', index=df.index, columns='gtex', aggfunc='first')
df.unstack()
df["New"] = df.groupby(["Note"]).cumcount()+1
df = df.set_index(["Note", "New"]).unstack().dropna(axis=1).sort_index()

#to_remove = ["153,full_data/ChIAPET_hg38_HT1376_cohesin_HT1376_hiseq_pairs.hic/", "153,full_data/ChIAPET_hg38_KU19_cohesin_KU19_hiseq_pairs.hic/", "43,full_data/ChIAPET_hg38_DU145_cohesin_DU145_hiseq_pairs.hic/", "53,full_data/ChIAPET_hg38_Jurkat_cohesin_Jurkat_hiseq_pairs.hic/"]
to_remove_includes = ["HiC_"]

#df = df.drop(index=to_remove)
for tr in to_remove_includes:
    df = df[~df.index.str.contains(tr)]

df.to_csv("final_results.csv")

cell_lines_mapping = {
22: "Cells_EBV-transformed_lymphocytes",
35: "Liver",
36: "Lung",
43: "Prostate",
50: "Thyroid",
53: "Whole_Blood",
69: "hepatocyte",
70: "neural progenitor cell",
77: "K562",
78: "HepG2",
79: "MCF-7",
93: "induced pluripotent stem cell",
114: "fibroblast of dermis",
115: "mammary epithelial cell",
122: "dermis microvascular lymphatic vessel endothelial cell",
144: "camera-type eye",
153: "urinary bladder",
163: "H1_Cell_Line",
166: "H1_Derived_Mesenchymal_Stem_Cells",
182: "Penis_Foreskin_Fibroblast_Primary_Cells_skin02",
210: "GM12878",
}
all_data = {}

for index, row in df["spearman"].iterrows():
    celline_id = int(index.split(",")[0])
    if not(celline_id in all_data):
        all_data[celline_id] = {}
    dataset_name = index.split("/")[1].split("/")[0]
    factor = "Baseline"
    if("cohesin" in dataset_name):
        factor = "Cohesin"
    elif("CTCF" in dataset_name):
        factor = "CTCF"
    elif("RNAPOL2" in dataset_name):
        factor = "RNAPOL2"
    dataset_name = dataset_name.removeprefix("ChIAPET_hg38_")
    dataset_name = dataset_name.removeprefix("ChIA_PIPE_")
    
    if(factor == "Baseline"):
        dataset_name = "Baseline"
    else:   
        seps = ['_cohesin', '_RNAPOL2', '_CTCF']
        for sep in seps:
            if(sep in dataset_name):
                dataset_name = factor + " ChIA-PET " + dataset_name.split(sep, 1)[0]
    all_data[celline_id][dataset_name] = row.to_list()

rows = 8
columns = 4

fig, axs = plt.subplots(rows, columns,  figsize=(21, 30))

i = 0
j = 0

def avg(arr):
    return sum(arr)/len(arr)

def create_single_chart(ax, baseline, dataset, dataset_name, title):
    df = pd.DataFrame({"Category": ["Baseline" for b in baseline] + [dataset_name for b in dataset], "Spearman correlation score": baseline + dataset})
    
    plt.sca(ax)
    ax = sns.stripplot(x='Category', y='Spearman correlation score', data=df, palette=my_palette, linewidth=1.5)
    ax = sns.boxplot(x='Category', y='Spearman correlation score',  palette=my_palette, data=df, showfliers=False)
    ax.set_ylabel("Spearman correlation score", fontsize=5)
    ax.set_xlabel(" ", fontsize=5)
    ax.set_title(title, fontsize=8, y = 1.1)
    _, xlabels = plt.xticks()
    ax.set_yticklabels(ax.get_yticks(), size = 5)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%0.3f'))
    ax.set_xticklabels(xlabels, size = 5)    
    return add_stat_annotation(ax, data=df, x='Category', y='Spearman correlation score',
                        box_pairs=[df["Category"].unique().tolist()], test='t-test_welch', text_format='star', loc='outside', verbose=2)

metrics = {}

for celline_id, datasets in all_data.items():
    baseline = datasets["Baseline"]
    for dataset_name, dataset in datasets.items():
        if(dataset_name == "Baseline"):
            continue
        if(avg(dataset) >= avg(baseline)):
            my_palette = ["#626567", "#2ECC71"]
        else:
            my_palette = ["#626567", "#E74C3C"]
        stat_metrics = create_single_chart(axs[j,i], baseline, dataset, dataset_name, cell_lines_mapping[celline_id])
        metrics[dataset_name] = {"pvalue": stat_metrics[1][0].pval, "dataset_corr": avg(dataset), "baseline_corr": avg(baseline), "diff_corr": avg(dataset)-avg(baseline)}
        i += 1
        if(i > columns-1):
            i = 0
            j += 1
            
plt.tight_layout()
plt.savefig("full_results.png", dpi=400)

factors = {"Cohesin" : {"baseline": [], "dataset": [] }, "CTCF" : {"baseline": [], "dataset": [] }, "RNAPOL2":  {"baseline": [], "dataset": [] }, "All factors":  {"baseline": [], "dataset": [] }}

fig, axs = plt.subplot_mosaic("ABC;DDD;DDD", figsize=(10, 10))
charts = ["A", "B", "C", "D"]
i = 0
for factor in factors.keys():
    for celline_id, datasets in all_data.items():
        baseline = datasets["Baseline"]
        for dataset_name, dataset in datasets.items():
            if(dataset_name == "Baseline"):
                continue
            if(factor in dataset_name):
                factors[factor]["baseline"] += baseline
                factors[factor]["dataset"] += dataset
                
                factors["All factors"]["baseline"] += baseline
                factors["All factors"]["dataset"] += dataset
                
    stat_metrics = create_single_chart(axs[charts[i]], factors[factor]["baseline"], factors[factor]["dataset"], factor, "Baseline vs " + factor)
    metrics[factor] = {"pvalue": stat_metrics[1][0].pval, "dataset_corr": avg(factors[factor]["dataset"]), "baseline_corr": avg(factors[factor]["baseline"]), "diff_corr": avg(factors[factor]["dataset"])-avg(factors[factor]["baseline"])}
    i+= 1
plt.tight_layout()
plt.savefig("factor_results.png", dpi=400)
pd.DataFrame(metrics).transpose().to_csv("full_results_metrics.csv")