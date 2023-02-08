import pandas as pd

df = pd.read_csv("5kbp_results.csv", index_col=None)[["Note", "spearman"]]

#df = df.pivot_table(values='spearman', index=df.index, columns='gtex', aggfunc='first')
df.unstack()
df["New"] = df.groupby(["Note"]).cumcount()+1
df = df.set_index(["Note", "New"]).unstack().dropna(axis=1).sort_index()

df.to_csv("5kbp_results_parsed.csv")

pass