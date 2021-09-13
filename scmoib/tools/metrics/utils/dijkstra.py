from multiprocessing import Process, Queue
from collections import defaultdict
from heapq import *
import numpy as np


def __dijkstra(edges, f, t):
    """Finds the shortest path between two nodes by using binary heap version of Dijkstra algorithm.
    
    Parameters
    ----------
    edges: list
        List of graph's edges.
    f: str
        Start node.
    t: str
        Target node.
    """
    g = defaultdict(list)
    for l, r, c in edges:
        g[l].append((c, r))

    q, seen, mins = [(0, f, ())], set(), {f: 0}
    while q:
        (cost, v1, path) = heappop(q)
        if v1 not in seen:
            seen.add(v1)
            path += (v1,)
            if v1 == t:
                return cost, path

            for c, v2 in g.get(v1, ()):
                if v2 in seen:
                    continue
                prev = mins.get(v2, None)
                next = cost + c
                if prev is None or next < prev:
                    mins[v2] = next
                    heappush(q, (next, v2, path))

    return float("inf"), None


def __adj_matr_to_edges(adata):
    """Converts adjacency matrix to edges list
    
    Parameters
    ----------
    adj_matr: np.array
        Graph's adjacency matrix
    
    Returns
    -------
    edges: list
        List of graph's edges
    
    """
    adj_matr = adata.obsp['connectivities'].A
    inds = np.where(adj_matr > 0)
    bc_list = list(adata.obs.index)
    edges = []
    for i, j in zip(inds[0], inds[1]):
        edges.append((bc_list[i], bc_list[j], adj_matr[i, j]))
    return edges


def __batch_generator(list1, list2, num_processes):
    """
    None
    """
    batch_size = len(list1) // num_processes
    for i in range(num_processes):
        ind1, ind2 = batch_size * i, batch_size * (i + 1)
        if i == num_processes - 1:
            ind2 = len(list1)
        yield list1[ind1:ind2], list2[ind1:ind2]


def __calculate_paths(edges, list1, list2, queue):
    """
    Function for multiprocessing step
    """
    res = []
    for i, j in zip(list1, list2):
        tmp_res = __dijkstra(edges, i, j)
        res.append((tmp_res[0], tmp_res[1]))
    queue.put(res)


def run_dijkstra(adata, bc_list1, bc_list2, n_jobs=None):
    """
    Calculate shortest paths for all pairs of nodes.
    
    Parameters
    ----------
    adata: AnnData object
    
    bc_list1: list
    
    bc_list2: list
    
    n_jobs: None or int, optional
        The number of jobs to use for the computation. None means 1 
    """
    if not n_jobs or n_jobs <= 0:
        n_jobs = 1
    edges = __adj_matr_to_edges(adata)
    processes = []
    q = Queue()
    batches = __batch_generator(bc_list1, bc_list2, n_jobs)
    results = []

    for i in range(n_jobs):
        tmp_batch = next(batches)
        p = Process(target=__calculate_paths, args=(edges, tmp_batch[0], tmp_batch[1], q))
        p.start()
        processes.append(p)

    for _ in range(len(processes)):
        results += q.get()

    for p in processes:
        p.join()

    return results
