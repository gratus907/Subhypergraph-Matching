import sys
from subprocess import Popen, PIPE
import pandas as pd
import platform
from experiment_info import *
import numpy as np
import tqdm 


class Result:
    def __init__(self):
        self.query_name = ""
        self.elapsed_time = 0.0
        self.num_results = 0
        self.num_recursion_calls = 0

    def __str__(self):
        return f"{self.query_name}: {self.num_results} results in {self.elapsed_time} ms (Rec = {self.num_recursion_calls})"

    def to_dict(self):
        return {
            "query_name": self.query_name,
            "num_results": self.num_results,
            "num_recursion_calls": self.num_recursion_calls,
            "elapsed_time": self.elapsed_time,
        }


def load_results(path : str):
    L = open(path).readlines()
    answers = {}
    for l in L:
        answers[l.split()[0]] = int(l.split()[1])
    return answers

def generate_args(binary, *params):
    arguments = [binary]
    arguments.extend(list(params))
    return arguments


def execute_binary(args):
    cmd = ' '.join(args)
    print(cmd)
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=sys.stderr)
    (std_output, std_error) = process.communicate()
    process.wait()
    rc = process.returncode
    return rc, std_output, std_error


def run(dataset):
    dataset_name = dataset.split('_')[0]
    queries_path = f"SubgraphMatching/{dataset_name}/"+(f"{dataset}_ans.txt")
    answers = load_results(queries_path)
    num_queries = len(answers)
    cmd = f"../build/SubgraphMatching -d {dataset_name} -q {queries_path}" 
    print(cmd)
    results : list[dict] = []
    total_time = 0.0
    r = Result()
    progress = tqdm.tqdm(total=num_queries)
    with Popen(cmd, shell=True, stdout=PIPE, stderr=sys.stderr, universal_newlines=True,bufsize=1) as process:
        for line in iter(process.stdout.readline, ''):
            if "Start" in line:
                r.query_name = line.split(" ")[-3].split('/')[-1].strip()
            if "Out:" in line:
                r.num_results = int(line.split(":")[1].strip())
            if "Rec:" in line:
                r.num_recursion_calls = int(line.split(":")[1].strip())
            if "Total:" in line:
                total_time = float(line.split(":")[1].strip())
            if "Time:" in line:
                r.elapsed_time = float(line.split(":")[1].strip())
                results.append(r.to_dict())
                r = Result()
                progress.update(1)
    process.wait()
    print(f"\nTotal time: {total_time} ms")
    for result in results:
        qname = result['query_name']
        if answers[qname] != result['num_results']:
            raise AssertionError(f"Ans for {qname} should be {answers[qname]} but is {result['num_results']}")
    print(f"Correctness Asserion OK")
    return pd.DataFrame(results)
    

do_not_save = False
exp_record = pd.read_csv("SubgraphMatching/test-results.csv")
gitver = get_current_repository_ver()[-8:]
currentdate = get_experiment_date()
dataset = sys.argv[1]
if dataset is None:
    dataset = 'yeast'
if (dataset, currentdate, gitver) in zip(exp_record['dataset'], exp_record['date'], exp_record['gitver']):
    print(f"Existing experiment metadata: {(dataset, currentdate, gitver)}")
    opt = input("Run without saving? ").lower()
    if 'yes' in opt: do_not_save = True  
    else: raise AssertionError(f"Existing experiment metadata: {(dataset, currentdate, gitver)}")
print(get_system_info())
print((dataset, currentdate, gitver))
results = run(dataset)

result_summary = (f"{dataset},{currentdate},{gitver},{results['elapsed_time'].mean():.3f},{results['num_recursion_calls'].mean():.3f}\n")
if do_not_save:
    print(result_summary)
else:
    result_path = open("SubgraphMatching/test-results.csv", 'a')
    result_path.write(result_summary)
