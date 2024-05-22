import sigProfilerPlotting as sigPlt


matrix_path = './test.csv'
output_path = './image2.pdf'
project = 'test'

sigPlt.plotSBS(matrix_path, output_path, project, '96', percentage=True)


from pymutspec.annotation import rev_comp
d = pd.read_csv("../data/MutSpecALLvert.csv")
d.columns = ['MutationType', 'Vert']
# d['Vert'] = (d['Vert'] * 10000).astype(int)
first_muts = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G',]

d1 = d[d.MutationType.str.slice(2, 5).isin(first_muts)]

d2 = d[~d.MutationType.str.slice(2, 5).isin(first_muts)]
d2['MutationType'] = d2.MutationType.apply(rev_comp)

d1.to_csv('../test_vert1_96.txt', sep='\t', index=False)
d2.to_csv('../test_vert2_96.txt', sep='\t', index=False)
import sigProfilerPlotting as sigPlt
matrix_path = '../test_vert1_96.txt'
output_path = '../'
project = 'vert_1'
sigPlt.plotSBS(matrix_path, output_path, project, '96', percentage=True)
import sigProfilerPlotting as sigPlt
matrix_path = '../test_vert2_96.txt'
output_path = '../'
project = 'vert_2'
sigPlt.plotSBS(matrix_path, output_path, project, '96', percentage=True)