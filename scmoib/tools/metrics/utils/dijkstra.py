from multiprocessing import Process, Queue
from collections import defaultdict
from heapq import *
import numpy as np
import warnings
from anndata import AnnData
from typing import Union, Tuple, List


def __dijkstra(edges: List[Tuple[str, str, float]], f: str, t: str):
    """Finds the shortest path between two nodes by using binary heap version of Dijkstra algorithm.
    
    Parameters
    ----------
    edges
        List of graph's edges.
    f
        Start node.
    t
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


def __adj_matr_to_edges(adata: AnnData) -> List[Tuple[str, str, float]]:
    """Converts adjacency matrix into the list of edges.
    
    Parameters
    ----------
    adata
        Anndata object containing connectivity matrix
    
    Returns
    -------
    edges
        List of graph's edges with their weights
    """
    adj_matr = adata.obsp['connectivities'].A
    inds = np.where(adj_matr > 0)
    bc_list = list(adata.obs.index)
    edges = []
    for i, j in zip(inds[0], inds[1]):
        edges.append((bc_list[i], bc_list[j], adj_matr[i, j]))
    return edges


def __batch_generator(list1: List[str], list2: List[str], num_processes: int):
    """Batch generator for multiprocessing step.

    Parameters
    ----------
    list1
        RNA matching barcodes
    list2
        ATAC matching barcodes
    num_processes
        The number of jobs to use for the computation.
    """
    batch_size = len(list1) // num_processes
    for i in range(num_processes):
        ind1, ind2 = batch_size * i, batch_size * (i + 1)
        if i == num_processes - 1:
            ind2 = len(list1)
        yield list1[ind1:ind2], list2[ind1:ind2]


def __calculate_paths(edges: List[Tuple[str, str, float]], list1: List[str], list2: List[str], queue: Queue):
    """Wrapper for multiprocessing step.

    Parameters
    ----------
    edges
        List of graph's edges with their weights.
    list1
        RNA matching barcodes.
    list2
        ATAC matching barcodes.
    queue
        Queue for storing results.
    """
    res = []
    for i, j in zip(list1, list2):
        tmp_res = __dijkstra(edges, i, j)
        res.append((tmp_res[0], tmp_res[1]))
    queue.put(res)


def run_dijkstra(
        adata: AnnData,
        bc_list1: List[str],
        bc_list2: List[str],
        n_jobs: Union[int, None] = None
) -> List[Tuple[float, Tuple[str]]]:
    """
    Calculate shortest paths for all pairs of barcodes.
    
    Parameters
    ----------
    adata
        Annotated data matrix
    bc_list1
        RNA matching barcodes
    bc_list2
        ATAC matching barcodes
    n_jobs
        The number of jobs to use for the computation. None means 1.

    Returns
    -------
    results
        List of shortest paths with their weights
    """
    warnings.warn(
        UserWarning(
            'Warning! You are running the high load graph method.\n'
            'Test run information:\n'
            'CPU: 8 cores\n'
            'Wall time: 6 minutes\n'
            'RAM: 4 GB(?)\n'
            'Number of barcode pairs: 11k'
        )
    )

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
