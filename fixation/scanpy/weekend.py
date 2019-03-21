### conda install seaborn scikit-learn statsmodels numba pytables
### conda install -c conda-forge python-igraph louvain
### pip install scanpy

import pickle,numpy
import matplotlib,matplotlib.pyplot
import scanpy
scanpy.settings.verbosity=5

# 1. read data
print('reading data...')
adata=scanpy.read_10x_mtx('/Volumes/omics4tb2/alomana/projects/i18/results/both_aggregated2/outs/filtered_feature_bc_matrix',var_names='gene_symbols',cache=True)
adata.var_names_make_unique() 
print(adata)
print()

# 2. filter data
print('filter data...')
scanpy.pp.filter_cells(adata,min_genes=200)
scanpy.pp.filter_genes(adata,min_cells=3)
print()

# 2.1. retrieve mouse cells only
print('retrieve mouse cells only...')
cellIDs=adata.obs_names.tolist()
geneNames=adata.var_names.tolist()

jarFile='/Volumes/omics4tb2/alomana/projects/i18/src/cs/scanpy/species.cellIDs.run.006.pickle'
f=open(jarFile,'rb')
[mouseCellIDs,humanCellIDs,chimericCellIDs]=pickle.load(f)
f.close()

print(len(cellIDs),len(geneNames))
print('mouse',len(mouseCellIDs))
print('human',len(humanCellIDs))
print('chimeric',len(chimericCellIDs))

# slice in mouse cells
print('before slicing mouse cells...')
print(adata)
print()
adata=adata[mouseCellIDs,:]
print('after')
print(adata)

# slice in mouse genes
mouseGenes=[element for element in geneNames if element[:4] == 'mm10']
print()
print('before slicing mouse genes...')
print(adata)
print()
adata=adata[:,mouseGenes]
print('after')
print(adata)
print()

# 3. QC
print('QC...')
mitoGenes=[element for element in adata.var_names if element[:8] == 'mm10_mt-']
print(mitoGenes)

adata.obs['percent_mito']=numpy.sum(adata[:,mitoGenes].X,axis=1).A1/numpy.sum(adata.X,axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
scanpy.pl.violin(adata,['n_genes','n_counts','percent_mito'],jitter=0.4,multi_panel=True,save='plot.QC.0.pdf',show=False)


print(adata)
scanpy.pl.scatter(adata, x='n_counts', y='percent_mito',save='plot.QC.1.pdf',show=False)
scanpy.pl.scatter(adata, x='n_counts', y='n_genes',save='plot.QC.2.pdf',show=False)


adata = adata[adata.obs['n_counts'] < 30000, :]
adata = adata[adata.obs['percent_mito'] < 0.05, :]
print(adata)
scanpy.pl.scatter(adata, x='n_counts', y='percent_mito',save='plot.QC.3.pdf',show=False)
scanpy.pl.scatter(adata, x='n_counts', y='n_genes',save='plot.QC.4.pdf',show=False)
print()

# 4. normalization and log transform
print('normalization...')
scanpy.pp.normalize_per_cell(adata, counts_per_cell_after=30e3)
scanpy.pp.log1p(adata)
print()

# 5. select highly variable genes
print('selecting HVGs...')
adata.raw = adata

scanpy.pp.highly_variable_genes(adata,min_mean=0.0125,max_mean=3,min_disp=0.5)
scanpy.pl.highly_variable_genes(adata,save='plot.HVG.pdf',show=False)
adata = adata[:,adata.var['highly_variable']]
print(adata)
print()

# 6. regress out effects and scale
print('regress and scale...')
scanpy.pp.regress_out(adata,['n_counts', 'percent_mito'])
scanpy.pp.scale(adata,max_value=10)
print()

# 7. associate treatment to cell ID
print('associate treatment to cellID...')
#f=open('/Volumes/omics4tb2/alomana/projects/i18/src/cs/scanpy/fixedLabels.pickle','rb')
#cellIDsFixed=pickle.load(f)
#f.close()

#f=open('/Volumes/omics4tb2/alomana/projects/i18/src/cs/scanpy/freshLabels.pickle','rb')
#cellIDsFresh=pickle.load(f)
#f.close() 


cellIDs=adata.obs_names.tolist()
print(len(cellIDs))

cellConditions=[]
for cellID in cellIDs:
    if '-1' in cellID: 
        cellConditions.append('fixed')
    elif '-2' in cellID:
        cellConditions.append('fresh')
    else:
        raise ValueError('cellID not recognized')

adata.obs['treatment']=cellConditions
    
print(cellConditions,len(cellConditions))
print()

# 8. Visualization
print('visualization...')

# 8.1. PCA
print('\t PCA...')
scanpy.tl.pca(adata, svd_solver='arpack')
scanpy.pl.pca(adata,color='treatment',palette=['red','blue'],alpha=0.5,save='plot.PCA.pdf',show=False)
print()

# 8.2. UMAP
print('\t UMAP...')
scanpy.pp.neighbors(adata,n_neighbors=10,n_pcs=30)
scanpy.tl.umap(adata)
scanpy.pl.umap(adata,color='treatment',alpha=0.5,save='plot.UMAP.pdf',show=False)
print()

print('\t Louvain...')
scanpy.tl.louvain(adata)
scanpy.pl.umap(adata, color=['louvain'],alpha=0.5,save='plot.UMAP.Louvain.pdf',show=False)

print('DONE.')
