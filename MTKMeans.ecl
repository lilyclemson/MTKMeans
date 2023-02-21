IMPORT ML_Core as core;
IMPORT core.Types as types;
IMPORT Datasets as ds;
IMPORT Python3 as Python;
IMPORT Std.system.Thorlib;

NF := types.numericField;
rawds :=  ds.ds0;
rawds1 := project(rawds, TRANSFORM({INTEGER id, RECORDOF(LEFT)}, SELF.id := COUNTER, SELF := LEFT));
Core.ToField(rawds1, NFds);
// output(NFds);


numCol := [];
n := COUNT(rawds);
col := MAX(NFds, number);
// OUTPUT(col, NAMED('col'));

nDs := NFds(number in numCol);
dDs := NFds(number not in numCol);
// OUTPUT(nDs, NAMED('nDS'));
// OUTPUT(dDs, NAMED('dDs'));

emptyDS := DATASET([], NF);

rmLayout := RECORD
    DECIMAL a;
    DECIMAL b;
END;

strLayout := RECORD
    STRING text
END;
colPairP := RECORD
    DECIMAL col0;
    DECIMAL col1;
    DECIMAL total;
    DECIMAL col1Cnt;
    REAL p;
END;

distsLayout := RECORD
    STRING  col;
    DECIMAL valueI;
    DECIMAL valueJ;
    DECIMAL dist;
    INTEGER n;
END;

STREAMED DATASET(distsLayout ) distances(STREAMED DATASET(NF) nDS = emptyDS, STREAMED DATASET(NF) dDS = emptyDS) := EMBED(Python: activity)
#STREAMED DATASET(strLayout) distances(STREAMED DATASET(NF) nDS = emptyDS, STREAMED DATASET(NF) dDS = emptyDS) := EMBED(Python: activity)
    import numpy as np
    import pandas as pd
    from itertools import combinations, permutations

    # set random seed
    np.random.seed(10)

    # read inputs
    dict = {}
    for i in dDS:
        wi, id, number, data = i
        if id in dict:
            dict[id].append(data)
        else:
            dict[id] = [data]
    dx = np.array(list(dict.values()))
    r,n = dx.shape
    

    ddf = pd.DataFrame(dx)
    ddf.columns = [str(i) for i in range(0, ddf.shape[1])]
    # All possible pairs of features
    colPairNames = permutations(ddf.columns, 2)
    
    # conditional probabilities of each attribute value in the first feature based on second feature 
    dictPij = {}
    for colPairName in list(colPairNames): # colPairName is tuple
        c0 = colPairName[0]
        c1 = colPairName[1]
        colPair = ddf[[c0, c1]]
        valuePair0 = colPair.groupby([c0, c1])[c1].count().reset_index(name='count')
        valuePair1 = colPair.groupby([c0]).agg(np.size)
        join = pd.merge(valuePair0, valuePair1, on = c0, how = 'left')
        join['p'] = join['count'] / join[c1 + '_y']
        '''    
        #yield(str(type(join)))
        npColPair = np.array(join)
        for i in range(len(npColPair)):
            yield((npColPair[i,0], npColPair[i,1], npColPair[i,2],npColPair[i,3],npColPair[i,4]))
            #yield(list(npColPair[i]))
    
        '''

        #yield(colPairName)
        dictPijKey = (c0, c1)
        col = {}
        for row in join.values.tolist():
            i = int(row[0])
            j = int(row[1]) # col J
                        
            if j in col:
                jCol = col[j]
                Pij = [i, row[4]]
                jCol.append(Pij)
                col[j] = jCol
            else:
                Pij = [[i, row[4]]]
                col[j] = Pij
        
        dictPij[dictPijKey] = col
        
        
        # construct shell dataframe for final result
        dists = {}
        for colName in ddf.columns:
            valuePairName = {}
            comb = combinations(ddf[colName].unique(), 2)
            for c in comb:
                valuePairName[c] = 0     # valuePair
            dists[i] = valuePairName
        

        
        for distsK, distsI in dists.items(): # key is the attribute name, item is the unique ValuePaires of the key
            #yield(str(type(key)))
            outKey = distsK
            subDistances = []
            for dictPijK, dictPijI in dictPij.items():
                if int(dictPijK[0]) == int( outKey):
                    subDistances.append(dictPijK)
            
            for valuePair, distance in distsI.items():
                value0 = valuePair[0]
                value1 = valuePair[1]
                values = []
                for subDistance in subDistances:
                    subColPair = dictPij[subDistance]
                    for subColPairK, subColPairI in subColPair.items():
                        larger = 0                            
                        for p in subColPairI :
                            if p[0] == value0 or p[0] == value1:                                    
                                if p[1] > larger:
                                    larger = p[1]
                        values.append(larger)
                d = float(sum(values)-(n-1))/(n-1)
                yield (str(outKey), float(value0), float(value1), d, n)  # (attribute, one attributeValue, another attributeValue, distance between attibute values)
        '''

        '''



ENDEMBED;


o := distances(nDs, dDs);
OUTPUT(o);



STREAMED DATASET(distsLayout) norm(STREAMED DATASET(NF) nDS = emptyDS, STREAMED DATASET(NF) dDS = emptyDS, UNSIGNED4 level = 2) := EMBED(Python: activity)
#STREAMED DATASET(strLayout) norm(STREAMED DATASET(NF) nDS = emptyDS, STREAMED DATASET(NF) dDS = emptyDS, UNSIGNED4 level = 2) := EMBED(Python: activity)
    import numpy as np
    
    # normalize Numeric attributes
    # read inputs
    dict = {}
    for i in nDS:
        wi, id, number, data = i
        if id in dict:
            dict[id].append(data)
        else:
            dict[id] = [data]
    nd = np.array(list(dict.values()))



    # min-max 
    nd = nd.astype(float)
    nd = (nd - nd.min(axis=0))/(x.max(axis=0) - x.min(axis=0))




ENDEMBED;

STREAMED DATASET(distsLayout) discretize(STREAMED DATASET(NF) dDS = emptyDS, UNSIGNED4 level = 2) := EMBED(Python: activity)
#STREAMED DATASET(strLayout) discretize(STREAMED DATASET(NF) dDS = emptyDS, UNSIGNED4 level = 2) := EMBED(Python: activity)
    import numpy as np

    # label discrete Data
    # read inputs
    dict = {}
    for i in dDS:
        wi, id, number, data = i
        if id in dict:
            dict[id].append(data)
        else:
            dict[id] = [data]
    dd = np.array(list(dict.values()))


    # digitize
    bins = [(1/level) * i for i in range(1, level)]
    dd = np.digitize(dd, bins = bins, right = True)
    



ENDEMBED;


initParams := RECORD
    UNSIGNED4 nodeID;
    UNSIGNED4 nNodes;
END;

STREAMED DATASET({UNSIGNED4 sessID}) init(STREAMED DATASET(initParams) initDat,
                                                        STRING wuid = WORKUNIT) := EMBED(Python: GLOBALSCOPE('mixedClustering'), 
                                                                                            PERSIST('global'), activity)

    import numpy as np
    global readInput, norm, discretize, initCenters

    def _readInput(ds, nCol = []):
        dictND = {}
        dictDD = {}
        for i in dDS:
            wi, id, number, data = i
            col = number -1
            if col in nCOl:
                if id in dictND:
                    dictND[id].append(data)
                else:
                    dictND[id] = [data]
            else:
                if id in dictDD:
                    dictDD[id].append(data)
                else:
                    dictDD[id] = [data]
        nd = np.array(list(dictND.values()))
        dd = np.array(list(dictDD.values()))
        
        return nd, dd
    

    def _norm(nd):

        # min-max 
        nd = nd.astype(float)
        result = (nd - nd.min(axis=0))/(x.max(axis=0) - x.min(axis=0))
        return result
    
    def _discretize(dd, level = 2): 

        # digitize
        bins = [(1/level) * i for i in range(1, level)]
        result = np.digitize(dd, bins = bins, right = True)
        return result

    def _initCenters(x, k, nCol=[]):
        row, col = x.shape
        

    readInput = _readInput
    norm = _norm
    discretize = _discretize
    for rec in initDat:
        nodeID, nNodes = rec 
        sessID = nNodes + int(wuid[1:9] + wuid[10:])
    return [(sessID, )]
    

ENDEMBED;

