
EXPORT Datasets := MODULE

EXPORT ds0 := DATASET([{1,1,1},{2,1,1}, {3,2,2}, {3,3,2},{4,3,3}, {4,2,3}], {INTEGER x, INTEGER y, INTEGER z});
SHARED l1 := RECORD
    STRING x;
    STRING y;
    REAL z;
END;
EXPORT ds1 := DATASET([{'A', 'C', 1.1},{'A', 'C', 2.9 },{'A', 'D', 3.1 },{ 'B', 'D', 4.9},{ 'B', 'C', 4.0},{'B', 'D', 3.2},{'A', ' D', 4.8}], l1);


END;