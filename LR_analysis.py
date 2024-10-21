# %%
import scipy.stats as sps
import pandas as pd
import numpy as np
import itertools as itt
from collections import Counter

data_file = "../data/malformation_exp_updated.ods"
sheet_name, wt_exp = "dand5", "14hpf_wt"
# sheet_name, wt_exp = "hearts", "2dpf_wt"
n_resamples = int(1e5)

df = pd.read_excel(data_file, sheet_name=sheet_name, index_col=0)
df_pval = df.copy()
df
# %%

categories = df.columns.to_list()
Nwt = df.loc[wt_exp, categories].to_numpy()


def statistic(N1, N2):
    P1 = N1 / N1.sum()
    P2 = N2 / N2.sum()
    return np.abs(P1 - P2).sum()


def resample(Ntot, n1):
    x = np.random.choice(range(sum(Ntot)), n1, replace=False)
    n_thr = 0
    N1 = np.zeros_like(Ntot)
    for i, n in enumerate(Ntot):
        N1[i] = ((n_thr <= x) & (x < n_thr + n)).sum()
        n_thr += n
    return N1, Ntot - N1


# %%
print(f"data file: {data_file}")
print(f"sheet name: {sheet_name}")
print(f"wild-type experiment: {wt_exp}")
print(f"categories: {categories}")
print(f"number of resamples: {n_resamples}")
print("------")

pval_dict = {}
for exp in df.index.to_list():
    print(f"processing {wt_exp} vs {exp}")

    Ns = df.loc[exp, categories].to_numpy()
    Ntot = Ns + Nwt
    s0 = statistic(Nwt, Ns)
    n_weird = 1
    for i in range(n_resamples):
        N1, N2 = resample(Ntot, Nwt.sum())
        assert np.all(N1 + N2 == Ntot)
        s = statistic(N1, N2)
        if s >= s0:
            n_weird += 1
    pval_dict[exp] = n_weird / (n_resamples + 1)
    print(f"pval <= {pval_dict[exp]:.2}")


# %%

df_pval[f"pval_{wt_exp}"] = df_pval.index.map(pval_dict)
df_pval.to_csv(f"../results/pvals/{sheet_name}_{wt_exp}_reshuffle.csv")

# %%