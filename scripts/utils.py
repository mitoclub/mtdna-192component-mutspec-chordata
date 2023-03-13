import numpy as np
import pandas as pd


def calc_mutspec_class(df: pd.DataFrame, group_col="Class"):
    ms_cls = df.groupby([group_col, 'Mut'])['RawMutSpec'].sum().reset_index()
    ms_cls["RawMutSpecSum"] = ms_cls.Class.map(ms_cls.groupby("Class").RawMutSpec.sum().to_dict())
    ms_cls['MutSpec'] = ms_cls.RawMutSpec / ms_cls.RawMutSpecSum
    ms_cls = ms_cls.drop(['RawMutSpec', 'RawMutSpecSum'], axis=1)
    return ms_cls


effect2sbs = {
    'MMR_deficiency': ['SBS6', 'SBS14', 'SBS15', 'SBS20', 'SBS21', 'SBS26', 'SBS44', 'DBS7', 'DBS10', 'ID7'], 
    'POL_deficiency': ['SBS10a', 'SBS10b', 'SBS10c', 'SBS10d', 'SBS28', 'DBS3'], 
    'HR_deficiency': ['SBS3', 'ID6'], 
    'BER_deficiency': ['SBS30', 'SBS36'], 
    'Chemotherapy': ['SBS11', 'SBS25', 'SBS31', 'SBS35', 'SBS86', 'SBS87', 'SBS90', 'DBS5'], 
    'Immunosuppressants': ['SBS32'], 
    'Treatment': ['SBS11', 'SBS25', 'SBS31', 'SBS32', 'SBS35', 'SBS86', 'SBS87', 'SBS90', 'DBS5'], 
    'APOBEC': ['SBS2', 'SBS13'], 
    'Tobacco': ['SBS4', 'SBS29', 'SBS92', 'DBS2', 'ID3'], 
    'UV': ['SBS7a', 'SBS7b', 'SBS7c', 'SBS7d', 'SBS38', 'DBS1', 'ID13'], 
    'AA': ['SBS22'], 
    'Colibactin': ['SBS88', 'ID18'], 
    'Artifact': ['SBS27', 'SBS43', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59', 'SBS60', 'SBS95'], 
    'Lymphoid': ['SBS9', 'SBS84', 'SBS85']
}
sbs2effect = {sbs: ef for ef,sbs_lst in effect2sbs.items() for sbs in sbs_lst}
sbs2effect["SBS18"] = "ROS"
