import os, sys

import os,sys,numpy
import sklearn,sklearn.decomposition,sklearn.manifold
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

quantDir='/Volumes/omics4tb2/alomana/projects/i18/results/bulk_rnaseq/kallisto.160/host/'
files=os.listdir(quantDir)

# 3. generating full expression matrix
print('generating expression matrix file...')

# 3.1. reading the genes
genes=[]
oneFile=quantDir+files[0]+'/abundance.tsv'
f=open(oneFile,'r')
next(f)
for line in f:
    vector=line.split()
    geneName=vector[0]
    genes.append(geneName)
f.close()


# 3.2. reading expression
expression={}
for element in files:
    if element not in expression.keys():
        expression[element]={}

    workingFile=quantDir+element+'/abundance.tsv'
    f=open(workingFile,'r')
    next(f)
    for line in f:
        vector=line.split()
        geneName=vector[0]
        abundance=float(vector[-1])
        expression[element][geneName]=abundance
    f.close()


# 3.3. writing expression matrix
expressionFile='expressionMatrix.kallisto.txt'
conditionNames=list(expression.keys())

x=[]

for i in range(len(files)):
    
    x.append([])

  
for i in range(len(genes)):
    for j in range(len(files)):
        value=expression[files[j]][genes[i]]

        x[j].append(value)
        
    

# 4. exploring the data
print('visualizing the data...')
original=numpy.array(x)

# 4.1. PCA of samples
print('running PCA...')
pcaMethod=sklearn.decomposition.PCA(n_components=5)
pcaObject=pcaMethod.fit(original)
new=pcaObject.transform(original)
explainedVar=pcaObject.explained_variance_ratio_
print('cumsum explained variance...')
print(numpy.cumsum(explainedVar))



for i in range(len(new)):
    matplotlib.pyplot.scatter(new[i,0],new[i,1],c='black',marker='o',s=60,edgecolors='black')

matplotlib.pyplot.xlabel('PCA 1 ({0:.2f} var)'.format(explainedVar[0]))
matplotlib.pyplot.ylabel('PCA 2 ({0:.2f} var)'.format(explainedVar[1]))
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figure.pca.png')
matplotlib.pyplot.clf()
print()

# 4.2. t-SNE of samples
print('running t-SNE...')
tSNE_Method=sklearn.manifold.TSNE(method='exact',verbose=1,init='pca')
tSNE_Object=tSNE_Method.fit(original)
new=tSNE_Object.fit_transform(original)

for i in range(len(new)):
    matplotlib.pyplot.scatter(new[i,0],new[i,1],c='black',marker='o',s=60,edgecolors='black')
matplotlib.pyplot.savefig('figure.tSNE.png')
matplotlib.pyplot.clf()
print()

print('... all done.')
