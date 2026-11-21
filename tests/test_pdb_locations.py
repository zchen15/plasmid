import numpy as np
import pandas as pd
from plasmid.pdb.locations import Location

def test_eq():
    loc0 = Location(range(0, 10))
    loc1 = Location(range(0, 10))
    assert loc0 == loc1    

def test_add():
    loc1 = Location(range(0, 10))
    loc2 = Location(range(10, 20))
    loc3 = loc2 - 5
    loc4 = loc2 + 5
    print('all locs', loc1, loc2, loc3, loc4)
    
    # add contiguous locations
    x = loc1 + loc2
    print('loc1 + loc2 = ', x)
    assert x == Location(range(0, 20))

    # add overlapping locations
    x = loc1 + loc3
    print('loc1 + loc3 = ', x)
    assert x == Location(range(0, 15))
    
    # add discontinous locations
    x = loc1 + loc4
    y = Location(range(0, 10)) + Location(range(15, 25))
    print('loc1 + loc4 = ', x)
    assert x == y

def test_sub():
    loc1 = Location(range(10, 50))
    print('loc1 = ', loc1)
    
    # sub overlapping locations
    x = loc1 - Location(range(20, 30))
    print('loc1 - Location(20,30) = ', x)
    assert x == Location(range(10, 20)) + Location(range(30, 50))
    
    x = loc1 - Location(range(40, 60))
    print('loc1 - Location(40, 60) = ', x)
    assert x == Location(range(10, 40))
    
    x = loc1 - Location(range(5, 15))
    print('loc1 - Location(5, 15) = ', x)
    assert x == Location(range(15, 50))

def test_df_sum():
    loc1 = Location(range(0, 10))
    loc2 = Location(range(10, 20))
    loc3 = Location(range(30, 35))
    
    df1 = pd.DataFrame([loc1, loc2-5, loc3], columns=["atom_id"])
    df1["chain_id"] = "A"
    print(df1)
    
    df2 = pd.DataFrame([loc1, loc2+5, loc3], columns=["atom_id"])
    df2["chain_id"] = "B"
    print(df2)
    
    df = pd.concat([df1, df2])
    x = df.groupby(["chain_id"]).agg({"atom_id":"sum"}).reset_index()
    print(x)
    assert x.values[0][1] == loc1 + (loc2 - 5) + loc3
    assert x.values[1][1] == loc1 + (loc2 + 5) + loc3

    y = x.groupby(["chain_id"]).agg({"atom_id":"sum"}).reset_index()
    print(y)
    assert len(y) == 2

