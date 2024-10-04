import pandas as pd
from pymutspec.annotation import calculate_mutspec, CodonAnnotation


coda = CodonAnnotation(2)

outdir = './data/for_asymmetry'

# Load preprocessed mutations dataset and site-specific reference annotation
mutations = pd.read_csv("./data/mutations.csv")
Ref = pd.read_csv("./data/ref_annot.csv")

major_arc_min_pos = 5800
major_arc_max_pos = 16000

# # v1
region_size = 5000
# v2
# region_size = 2500
print("Split mutations dataset based on genome site (only major arc):")
print(f"Low TSSS region: {major_arc_min_pos}-{major_arc_min_pos+region_size}\n"
      f"High TSSS region: {major_arc_max_pos-region_size}-{major_arc_max_pos}\n")

# Prepare 2 samples of mutations: LOW and HIGH
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

# Prepare tables of SYN observed mutations only
obs_low_tsss_syn = obs_low_tsss_all[obs_low_tsss_all.Label >= 1]
obs_high_tsss_syn = obs_high_tsss_all[obs_high_tsss_all.Label >= 1]

print("Number of mutations:")
print("Low TSSS, SYN+NONSYN:", obs_low_tsss_all.shape[0])
print("High TSSS, SYN+NONSYN:", obs_high_tsss_all.shape[0])
print("Low TSSS, SYN:", obs_low_tsss_syn.shape[0])
print("High TSSS, SYN:", obs_high_tsss_syn.shape[0], end="\n")


# Prepare expected mutation freqs for spectra normalization
# Exclude D-loop region due to different genome functional properties
cur_ref = Ref[(Ref.Type != "D-loop")].assign(AltNuc="ACGT")
cur_ref["AltNuc"] = cur_ref.AltNuc.apply(list)
cur_ref["AltCodon"] = cur_ref.apply(lambda x: coda.get_syn_codons(x.Codon, x.PosInCodon-1) if x.PosInCodon > -1 else [], axis=1)

# Split genome to LOW and HIGH
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

# Prepare expected mutations for LOW
exp_low_tsss_all = ref_low_tsss.explode("AltNuc")
exp_low_tsss_all = exp_low_tsss_all[exp_low_tsss_all.Nuc != exp_low_tsss_all.AltNuc]
exp_low_tsss_all["Sbs12"] = exp_low_tsss_all.Nuc + ">" + exp_low_tsss_all.AltNuc
exp_low_tsss_all["Sbs192"] = exp_low_tsss_all.Context.str.get(0) + "[" + exp_low_tsss_all["Sbs12"] + "]" + exp_low_tsss_all.Context.str.get(-1)

# Same for SYN expected mutations
exp_low_tsss_syn = ref_low_tsss.explode("AltCodon").dropna(subset="AltCodon")
exp_low_tsss_syn["Sbs12"] = exp_low_tsss_syn.Nuc + ">" + exp_low_tsss_syn.apply(lambda x: x.AltCodon[x.PosInCodon-1], axis=1)
exp_low_tsss_syn["Sbs192"] = exp_low_tsss_syn.Context.str.get(0) + "[" + exp_low_tsss_syn["Sbs12"] + "]" + exp_low_tsss_syn.Context.str.get(-1)

# Prepare expected mutations for HIGH
exp_high_tsss_all = ref_high_tsss.explode("AltNuc")
exp_high_tsss_all = exp_high_tsss_all[exp_high_tsss_all.Nuc != exp_high_tsss_all.AltNuc]
exp_high_tsss_all["Sbs12"] = exp_high_tsss_all.Nuc + ">" + exp_high_tsss_all.AltNuc
exp_high_tsss_all["Sbs192"] = exp_high_tsss_all.Context.str.get(0) + "[" + exp_high_tsss_all["Sbs12"] + "]" + exp_high_tsss_all.Context.str.get(-1)

# Same for SYN expected mutations
exp_high_tsss_syn = ref_high_tsss.explode("AltCodon").dropna(subset="AltCodon")
exp_high_tsss_syn["Sbs12"] = exp_high_tsss_syn.Nuc + ">" + exp_high_tsss_syn.apply(lambda x: x.AltCodon[x.PosInCodon-1], axis=1)
exp_high_tsss_syn["Sbs192"] = exp_high_tsss_syn.Context.str.get(0) + "[" + exp_high_tsss_syn["Sbs12"] + "]" + exp_high_tsss_syn.Context.str.get(-1)

# Calculate region- and muttype-specific freqs of expected mutations
exp_low_tsss_all_freqs = exp_low_tsss_all.Sbs192.value_counts().to_dict()
exp_low_tsss_syn_freqs = exp_low_tsss_syn.Sbs192.value_counts().to_dict()
exp_high_tsss_all_freqs = exp_high_tsss_all.Sbs192.value_counts().to_dict()
exp_high_tsss_syn_freqs = exp_high_tsss_syn.Sbs192.value_counts().to_dict()


# Calculate mutational spectra for 4 samples of mutations with specific expected values
ms_low_tsss_all = calculate_mutspec(obs_low_tsss_all, exp_low_tsss_all_freqs, use_context=True)
ms_low_tsss_syn = calculate_mutspec(obs_low_tsss_syn, exp_low_tsss_syn_freqs, use_context=True)
ms_high_tsss_all = calculate_mutspec(obs_high_tsss_all, exp_high_tsss_all_freqs, use_context=True)
ms_high_tsss_syn = calculate_mutspec(obs_high_tsss_syn, exp_high_tsss_syn_freqs, use_context=True)

# Save spectra
ms_low_tsss_all.sort_values("Mut").to_csv(f"{outdir}/ms_low_tsss_all.csv", index=False)
ms_low_tsss_syn.sort_values("Mut").to_csv(f"{outdir}/ms_low_tsss_syn.csv", index=False)
ms_high_tsss_all.sort_values("Mut").to_csv(f"{outdir}/ms_high_tsss_all.csv", index=False)
ms_high_tsss_syn.sort_values("Mut").to_csv(f"{outdir}/ms_high_tsss_syn.csv", index=False)


print("\nPrint number of NaN values in final spectra:")
print("low, all:", ms_low_tsss_all.MutSpec.isna().sum())
print("low, syn:", ms_low_tsss_syn.MutSpec.isna().sum())
print("high, all:", ms_high_tsss_all.MutSpec.isna().sum())
print("high, syn:", ms_high_tsss_syn.MutSpec.isna().sum())

print("\nPrint number of zero-values in final spectra:")
print("low, all:", (ms_low_tsss_all.MutSpec == 0).sum())
print("low, syn:", (ms_low_tsss_syn.MutSpec == 0).sum())
print("high, all:", (ms_high_tsss_all.MutSpec == 0).sum())
print("high, syn:", (ms_high_tsss_syn.MutSpec == 0).sum())
