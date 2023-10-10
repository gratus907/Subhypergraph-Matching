import numpy as np 
import random
from tqdm.auto import tqdm
class HyperGraph:
    hyperedges : list
    incidence_list : list 
    degree_sequence : list
    node_labels : list
    hyperedge_vertex_weights : list

class Walk:
    node_mappings : dict
    hyperedges : list
    node_labels : list
    num_vertices : int = 0
    

def to_string(W : Walk, H : HyperGraph):
    s = [f"{W.num_vertices} {len(W.hyperedges)}"]
    node_maps = list(W.node_mappings.items())
    s.append(" ".join(list(map(lambda x : str(H.node_labels[x]), W.node_mappings.keys()))))
    for h in W.hyperedges:
        s.append(",".join(sorted(list(map(lambda x : str(W.node_mappings[x]), h)))))
    return ("\n".join(s))


def read_hypergraph(dataset : str):
    H = HyperGraph()
    H.hyperedges = list(set(map(lambda x : frozenset(map(lambda y : int(y)-1, x.strip().split(','))), open(f"{path}/hyperedges-{dataset}.txt").readlines())))
    H.node_labels = list(map(lambda x : int(x.strip())-1, open(f"{path}/node-labels-{dataset}.txt").readlines()))
    num_vertices = 0
    for h in tqdm(H.hyperedges):
        for v in h:
            num_vertices = max(num_vertices, v + 1)
    H.incidence_list = [list() for _ in range(num_vertices)]
    H.degree_sequence = [0] * num_vertices
    for i, h in enumerate(tqdm(H.hyperedges)):
        for v in h:
            H.incidence_list[v].append(i)
            H.degree_sequence[v] += 1
    H.hyperedge_vertex_weights = [list() for _ in range(len(H.hyperedges))]
    for i in tqdm(range(len(H.hyperedges))):
        H.hyperedges[i] = np.array(list(H.hyperedges[i]))
        H.hyperedge_vertex_weights[i] = np.array(
            list(map(lambda x : H.degree_sequence[x]-1,
                H.hyperedges[i])), dtype=np.float32
        )
        H.hyperedge_vertex_weights[i] /= np.sum(H.hyperedge_vertex_weights[i])
    H.hyperedges = np.array(H.hyperedges, dtype=object)
    return H

def generate_walk(H : HyperGraph, num_edges : int):
    current_vertex = -1
    chosen_hyperedges = []
    while len(chosen_hyperedges) < num_edges:
        e_idx = -1
        while e_idx == -1 or (e_idx in chosen_hyperedges):
            e_idx = (
                random.randint(0, len(H.hyperedges)) if len(chosen_hyperedges) == 0 
                else np.random.choice(H.incidence_list[current_vertex])
            )
        chosen_hyperedges.append(e_idx)
        current_vertex = np.random.choice(H.hyperedges[e_idx], p=H.hyperedge_vertex_weights[e_idx])
    w = Walk()
    w.hyperedges = H.hyperedges[chosen_hyperedges]
    w.node_mappings = dict()
    for h in w.hyperedges:
        for v in h:
            if v not in w.node_mappings:
                w.node_mappings[v] = w.num_vertices
                w.num_vertices += 1
    return w


def write_walk(file_path : str, w : Walk, H : HyperGraph):
    print(to_string(w, H), file=open(file_path, 'w'))
            

dataset = "amazon-reviews"
path = f"../dataset/hypergraphs/{dataset}"
print("Reading hypergraph...")
H = read_hypergraph(dataset)
print("OK! Walk generation...")
for i in [3, 4, 5, 6, 7, 8]:
    for tc in tqdm(range(200)):
        w = generate_walk(H, i)
        write_walk(f"{path}/queries/query_{i}_{tc}.txt", w, H)