from subprocess import Popen, PIPE
import pandas as pd
import platform
from machine_info import get_system_info
import sys
import numpy as np
from natsort import natsorted
import datetime
import multiprocessing
import os 
from tqdm import tqdm

now = datetime.datetime.now()
TIME = f"{now.month}-{now.day}-{now.hour}-{now.minute}"


binary_path = "../build/SubhypergraphMatching"

def execute_binary(cmd):
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    (std_output, std_error) = process.communicate()
    process.wait()
    rc = process.returncode
    return rc, std_output, std_error

import re
def escape_ansi(line):
    ansi_escape = re.compile(r'(?:\x1B[@-_]|[\x80-\x9F])[0-?]*[ -/]*[@-~]')
    return ansi_escape.sub('', line)

def run_query(args):
    data, query = args
    cmd = f"{binary_path} -d {data} -q {query}"
    (rc, std_output, std_error) = execute_binary(cmd)
    if rc == 5:
        return {}
    result = str(std_error, encoding='utf-8')
    results = result.split('\n')
    exec_time, cs_v_init, cs_e_init, cs_v_after, cs_e_after, Vq, Eq, Aq, Vg, Eg, Ag, Vl, El = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    state = 0
    for line in results:
        line = escape_ansi(line)
        if "PatternGraph" in line:
            state += 1
            continue
        elif "DataGraph" in line:
            state += 1 
            continue
        elif "CS Statistics" in line:
            state += 1
            continue
        
        if state == 1:
            if 'V, E, TotalArity' in line:
                Vq, Eq, Aq = map(int, line.split('=')[1].strip().split(','))
            elif '#VLabel' in line:
                Vl, El = map(int, line.split('=')[1].strip().split(','))
        elif state == 2:
            if 'V, E, TotalArity' in line:
                Vg, Eg, Ag = map(int, line.split('=')[1].strip().split(','))
        elif state == 3:
            if 'vertices' in line:
                cs_v_init = int(line.split(':')[1].strip())
            elif 'hyperedges' in line:
                cs_e_init = int(line.split(':')[1].strip())
        elif state == 4:
            if 'vertices' in line:
                cs_v_after = int(line.split(':')[1].strip())
            elif 'hyperedges' in line:
                cs_e_after = int(line.split(':')[1].strip())
        if 'FilteringTime' in line:
            exec_time = float(line.split(':')[1].strip())

    return {
        'dataset': data,
        'query': query,
        'Vq': Vq,
        'Eq': Eq,
        'Aq': Aq,
        'Vg': Vg,
        'Eg': Eg,
        'Ag': Ag,
        'Vl': Vl,
        'El': El,
        'cs_v_init': cs_v_init,
        'cs_e_init': cs_e_init,
        'cs_v_after': cs_v_after,
        'cs_e_after': cs_e_after,
        'exec_time': exec_time,
    }


import itertools
if __name__ == '__main__':
    os.system("cd ../build && cmake .. && make")
    print(get_system_info())
    query_names = []
    datasets = ['contact-high-school']
    for sz in [3,4,5,6]:
        for idx in range(200):
            query_names.append(f"query_{sz}_{idx}")
    L = [(dataset, query) for dataset, query in itertools.product(datasets, query_names)]
    pool = multiprocessing.Pool(processes=4)
    results = list(tqdm(pool.imap_unordered(run_query, L), total=len(L)))
    df = pd.DataFrame(results)
    print(df)
    df.to_csv("results.csv", index=False)