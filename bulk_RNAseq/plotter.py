import os,sys,numpy
import sklearn,sklearn.decomposition,sklearn.manifold
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

###
### MAIN
###

# 0. user defined variables
quantDir='/Volumes/omics4tb2/alomana/projects/i18/results/bulk_rnaseq/kallisto.160/host/'
resultsDir='/Volumes/omics4tb2/alomana/projects/i18/results/bulk_rnaseq/expression.host.160/'

# 1. generating full expression matrix
print('generating expression matrix file...')

if os.path.exists(resultsDir) == False:
    os.mkdir(resultsDir)
    
items=os.listdir(quantDir)
files=[element for element in items]

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
        print(geneName,abundance,element)

        sys.exit()
    f.close()

# 3.3. writing expression matrix
expressionFile=resultsDir+'expressionMatrix.kallisto.txt'
conditionNames=list(expression.keys())

rbfConditions=[element for element in conditionNames if 'rbf' in element]
inverse=[element[::-1] for element in rbfConditions]
inverse.sort()
revertedRBF=[element[::-1] for element in inverse]

trnaConditions=[element for element in conditionNames if 'trna' in element]
inverse=[element[::-1] for element in trnaConditions]
inverse.sort()
revertedTRNA=[element[::-1] for element in inverse]

reverted=revertedTRNA+revertedRBF

x=[]
theEdgeColors=[]
theFaceColors=[]
theMarkers=[]

g=open(expressionFile,'w')

g.write('\t')
for i in range(len(reverted)):
    g.write('{}\t'.format(reverted[i]))
    
    x.append([])

    if 'tp.1' in reverted[i]:
        theEdgeColors.append('blue')
    elif 'tp.2' in reverted[i]:
        theEdgeColors.append('green')
    elif 'tp.3' in reverted[i]:
        theEdgeColors.append('orange')
    else:
        theEdgeColors.append('red')

    if 'trna' in reverted[i]:
        theFaceColors.append('w')
    else:
        theFaceColors.append(theEdgeColors[i])

    if 'rep.1' in reverted[i]:
        theMarkers.append('o')
    elif 'rep.2' in reverted[i]:
        theMarkers.append('s')
    else:
        theMarkers.append('^')
    
g.write('\n')

for i in range(len(genes)):
    g.write('{}\t'.format(synonyms[genes[i]]))
    for j in range(len(reverted)):
        value=expression[reverted[j]][genes[i]]
        g.write('{}\t'.format(value))

        x[j].append(value)
        
    g.write('\n')
    
g.close()

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
    matplotlib.pyplot.scatter(new[i,0],new[i,1],c=theFaceColors[i],marker=theMarkers[i],s=60,edgecolors=theEdgeColors[i])

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
    matplotlib.pyplot.scatter(new[i,0],new[i,1],c=theFaceColors[i],marker=theMarkers[i],s=60,edgecolors=theEdgeColors[i])
matplotlib.pyplot.savefig('figure.tSNE.png')
matplotlib.pyplot.clf()
print()

print('... all done.')
