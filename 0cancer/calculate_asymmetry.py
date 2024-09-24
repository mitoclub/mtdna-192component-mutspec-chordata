import pandas as pd
from pymutspec.annotation import calculate_mutspec, CodonAnnotation


coda = CodonAnnotation(2)

mutations = pd.read_csv('./data/mutations.csv')
Ref = pd.read_csv("./data/ref_annot.csv")

major_arc_min_pos = 5800
major_arc_max_pos = 16000

region_size = 5000
print(f"Low TSSS region: {major_arc_min_pos}-{major_arc_min_pos+region_size}\nHigh TSSS region: {major_arc_max_pos-region_size}-{major_arc_max_pos}\n")

# **Prepare OBS**
obs_low_tsss_all = mutations[
    (mutations.Type != "D-loop") &
    (mutations.Pos > major_arc_min_pos) &
    (mutations.Pos < major_arc_min_pos+region_size)
]
obs_high_tsss_all = mutations[
    (mutations.Type != "D-loop") &
    (mutations.Pos > major_arc_max_pos-region_size) &
    (mutations.Pos < major_arc_max_pos)
]

obs_low_tsss_syn = obs_low_tsss_all[obs_low_tsss_all.Label >= 1]
obs_high_tsss_syn = obs_high_tsss_all[obs_high_tsss_all.Label >= 1]

print(obs_low_tsss_all.shape, obs_high_tsss_all.shape)
print(obs_low_tsss_syn.shape, obs_high_tsss_syn.shape)
print(
    obs_low_tsss_all.merge(Ref[['Pos', 'TBSS']]).TBSS.mean().round(),
    obs_low_tsss_syn.merge(Ref[['Pos', 'TBSS']]).TBSS.mean().round(),
    obs_high_tsss_all.merge(Ref[['Pos', 'TBSS']]).TBSS.mean().round(),
    obs_high_tsss_syn.merge(Ref[['Pos', 'TBSS']]).TBSS.mean().round(),
)

# **Prepare EXP**
cur_ref = Ref[(Ref.Type != "D-loop")].assign(AltNuc="ACGT")
cur_ref["AltNuc"] = cur_ref.AltNuc.apply(list)
cur_ref["AltCodon"] = cur_ref.apply(lambda x: coda.get_syn_codons(x.Codon, x.PosInCodon-1) if x.PosInCodon > -1 else [], axis=1)

ref_low_tsss = cur_ref[
    (cur_ref.Pos > major_arc_min_pos) &
    (cur_ref.Pos < major_arc_min_pos+region_size) & 
    (cur_ref.Strand == 1)
]
ref_high_tsss = cur_ref[
    (cur_ref.Pos > major_arc_max_pos-region_size) &
    (cur_ref.Pos < major_arc_max_pos) &
    (cur_ref.Strand == 1)
]
exp_low_tsss_all = ref_low_tsss.explode("AltNuc")
exp_low_tsss_all = exp_low_tsss_all[exp_low_tsss_all.Nuc != exp_low_tsss_all.AltNuc]
exp_low_tsss_all["Sbs12"] = exp_low_tsss_all.Nuc + ">" + exp_low_tsss_all.AltNuc
exp_low_tsss_all["Sbs192"] = exp_low_tsss_all.Context.str.get(0) + "[" + exp_low_tsss_all["Sbs12"] + "]" + exp_low_tsss_all.Context.str.get(-1)

exp_low_tsss_syn = ref_low_tsss.explode("AltCodon").dropna(subset="AltCodon")
exp_low_tsss_syn["Sbs12"] = exp_low_tsss_syn.Nuc + ">" + exp_low_tsss_syn.apply(lambda x: x.AltCodon[x.PosInCodon-1], axis=1)
exp_low_tsss_syn["Sbs192"] = exp_low_tsss_syn.Context.str.get(0) + "[" + exp_low_tsss_syn["Sbs12"] + "]" + exp_low_tsss_syn.Context.str.get(-1)


exp_high_tsss_all = ref_high_tsss.explode("AltNuc")
exp_high_tsss_all = exp_high_tsss_all[exp_high_tsss_all.Nuc != exp_high_tsss_all.AltNuc]
exp_high_tsss_all["Sbs12"] = exp_high_tsss_all.Nuc + ">" + exp_high_tsss_all.AltNuc
exp_high_tsss_all["Sbs192"] = exp_high_tsss_all.Context.str.get(0) + "[" + exp_high_tsss_all["Sbs12"] + "]" + exp_high_tsss_all.Context.str.get(-1)

exp_high_tsss_syn = ref_high_tsss.explode("AltCodon").dropna(subset="AltCodon")
exp_high_tsss_syn["Sbs12"] = exp_high_tsss_syn.Nuc + ">" + exp_high_tsss_syn.apply(lambda x: x.AltCodon[x.PosInCodon-1], axis=1)
exp_high_tsss_syn["Sbs192"] = exp_high_tsss_syn.Context.str.get(0) + "[" + exp_high_tsss_syn["Sbs12"] + "]" + exp_high_tsss_syn.Context.str.get(-1)


exp_low_tsss_all_freqs = exp_low_tsss_all.Sbs192.value_counts().to_dict()
exp_low_tsss_syn_freqs = exp_low_tsss_syn.Sbs192.value_counts().to_dict()
exp_high_tsss_all_freqs = exp_high_tsss_all.Sbs192.value_counts().to_dict()
exp_high_tsss_syn_freqs = exp_high_tsss_syn.Sbs192.value_counts().to_dict()

print(exp_low_tsss_all.shape)
print(exp_low_tsss_syn.shape)
print(exp_high_tsss_all.shape)
print(exp_high_tsss_syn.shape)
print(
    exp_low_tsss_all.TBSS.mean().round(),
    exp_low_tsss_syn.TBSS.mean().round(),
    exp_high_tsss_all.TBSS.mean().round(),
    exp_high_tsss_syn.TBSS.mean().round(),
)
ms_low_tsss_all = calculate_mutspec(obs_low_tsss_all, exp_low_tsss_all_freqs, use_context=True)
ms_low_tsss_syn = calculate_mutspec(obs_low_tsss_syn, exp_low_tsss_syn_freqs, use_context=True)
ms_high_tsss_all = calculate_mutspec(obs_high_tsss_all, exp_high_tsss_all_freqs, use_context=True)
ms_high_tsss_syn = calculate_mutspec(obs_high_tsss_syn, exp_high_tsss_syn_freqs, use_context=True)
ms_low_tsss_all.to_csv("./asymmetry/ms_low_tsss_all.csv", index=False)
ms_low_tsss_syn.to_csv("./asymmetry/ms_low_tsss_syn.csv", index=False)
ms_high_tsss_all.to_csv("./asymmetry/ms_high_tsss_all.csv", index=False)
ms_high_tsss_syn.to_csv("./asymmetry/ms_high_tsss_syn.csv", index=False)
print(ms_low_tsss_all.MutSpec.isna().sum(),
    ms_low_tsss_syn.MutSpec.isna().sum(),
    ms_high_tsss_all.MutSpec.isna().sum(),
    ms_high_tsss_syn.MutSpec.isna().sum()
)
print((ms_low_tsss_all.MutSpec == 0).sum(),
    (ms_low_tsss_syn.MutSpec == 0).sum(),
    (ms_high_tsss_all.MutSpec == 0).sum(),
    (ms_high_tsss_syn.MutSpec == 0).sum()
)