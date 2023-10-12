import pandas as pd
import sys
import os
from pathlib import Path

"""
TODO
"""
def get_sequential_data(benchmarks, file_path, timeout):
    seq_data = pd.read_csv(file_path)
    tts = []
    
    for b in benchmarks:
        solve_time = 2*timeout
        res = list(seq_data[seq_data["benchmark"] == b]["result"])[0]
        tim = list(seq_data[seq_data["benchmark"] == b]["time_real"])[0]
        if res != "unknown" and tim <= timeout:
            solve_time = tim
        tts.append(solve_time)
    return tts


"""
TODO
"""
def build_scamble_dictionary(benchmarks, scram_logfile, scram_loc):
    log_colms=['benchmark', 'time_real'] 
    scram_log_data = pd.read_csv(scram_logfile, names=log_colms)

    scram_dict = {}

    for b in benchmarks:
        name = Path(b).stem
        csvfile = (f"{scram_loc}{name}_combined.csv")
        data = pd.read_csv(csvfile)
        scram_dict[b] = data

    return scram_dict, scram_log_data

"""
TODO
"""
def access_scram_data(j, sname, data, scram_logdata):
    # get run time
    res = list(data[data["benchmark"] == sname]["result"])[0]
    tim = list(data[data["benchmark"] == sname]["time_real"])[0]
        
    # get scrambling time, skip s0
    if j == 0:
        return res, tim
    
    scram_time = (
        list(scram_logdata[scram_logdata["benchmark"].str.contains(sname)]["time_real"])[0])

    tim += scram_time
    return res, tim

"""
TODO
"""
def get_scrambled_data(benchmarks, n, scramdata, scram_logdata, timeout, 
                       prev_results=[]):
    tts = []
    
    i = 0
    for b in benchmarks:
        name = Path(b).stem
        try:
            data = scramdata[b]
            
            if len(prev_results) == 0:
                min_time = 2*timeout
                for j in range(n):
                    sname = f"{name}_s{j}.smt2"
                    res, tim = access_scram_data(j, sname, data, 
                                                 scram_logdata)
                    if res != "unknown":
                        if tim < min_time and tim <= timeout:
                            min_time = tim
            else:
                min_time = prev_results[i]
                j = n-1
                sname = f"{name}_s{j}.smt2"
                res, tim = access_scram_data(j, sname, data, 
                                                 scram_logdata)
                if res != "unknown":
                    if tim < min_time  and tim <= timeout :
                        min_time = tim
                
            tts.append(min_time)
            i += 1
        except Exception as e:
            print("exception", e)
            i += 1
            continue
    return tts

"""
TODO
"""
def get_num_solved(lst):
    su = 0
    for l in lst:
        if l < 2400:
            su += 1
    return su

"""
TODO
"""
def get_unknown_time(df):
    tts = 0
    tts_partitions = []
    if "sat" in list(df.result):
        tts = min(list(df[df.result == "sat"].time_real))
        tts_partitions = list(df.time_real)
    elif "unknown" in list(df.result):
        # benchmark was unsolved, oh well
        tts = 2400
        tts_partitions = list(df.time_real) 
    # if no sat found and no unknown found, 
    # then all are unsat, hooray!
    else:
        tts = max(list(df.time_real))
        tts_partitions = list(df.time_real)
    return tts

### For x axis 
"""
TODO
"""
def forward(x):
    return x**(1/2)

"""
TODO
"""
def inverse(x):
    return x**2

