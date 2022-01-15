import gzip
import pickle
import numpy as np
from sklearn.neighbors import NearestNeighbors

from ordered_set import OrderedSet

import bz2
import pickle 
import _pickle as cPickle


def transform(vecs):
    maxnorm = max([np.linalg.norm(v) for v in vecs])
    new_vecs = []
    for v in vecs:
        new_vecs.append(np.insert(v, 0, np.sqrt(maxnorm**2 - np.linalg.norm(v)**2)))

    return new_vecs


def append_one_base(kmers):
    new_kmers = []
    for kmer in kmers:
        for b in ['A', 'C', 'G', 'T']:
            new_kmers.append(f'{kmer}{b}')

    return new_kmers


def make_kmers(k):
    kmers = [""]
    for i in range(k):
        kmers = append_one_base(kmers)

    return kmers


def make_kmer_to_pos_dict(kmers):
    kpdict = {}
    for i,k in enumerate(kmers):
        kpdict[k] = i

    return kpdict


def make_vector(bc, kpdict):
    k = len(next(iter(kpdict)))

    vec = [0] * 4**k
    for i in range(len(bc) - k + 1):
        vec[kpdict[bc[i:i+k]]] += 1

    return vec


def load_barcode_whitelist(filename):
    barcodes = OrderedSet()

    if filename is not None:
        f = gzip.open(filename, 'rb') if filename.endswith(".gz") else open(filename, 'rb')
        for l in f:
            bc = l.decode("utf-8").strip()

            barcodes.add(bc)

    # for i in range(5):
    #     q = make_vector(f'A{barcodes[i]}C', kpdict)
    #     q_trans = np.insert(q, 0, 0)

    #     start = time.time()
    #     distances, indices = nbrs.kneighbors(np.array([q_trans]))
    #     print("Elapsed time", time.time() - start)
    #     print("Min index", indices[0])

    #     for i in range(len(indices[0])):
    #         # print(f'{q} {vecs[indices[0][i]]}')
    #         print(f'{barcodes[i]} {barcodes[indices[0][i]]}')

    #     print("")

    # print("")

    # a = pa.msa_aligner()
    # res=a.msa(barcodes, out_cons=False, out_msa=False, incr_fn='') # perform multiple sequence alignment 
    # res.print_msa() # print row-column multiple sequence alignment in PIR format

    return barcodes


def prepare_database(barcodes, k=3, n_neighbors=10):
    kpdict = make_kmer_to_pos_dict(make_kmers(k))

    vecs = [make_vector(bc, kpdict) for bc in barcodes]
    vecs_trans = transform(vecs)

    X = np.array(vecs_trans)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='kd_tree').fit(X)

    return nbrs


def save_database(knn, knnfile):
    with bz2.BZ2File(knnfile, 'w') as f:
        cPickle.dump(knn, f)


def load_database(knnfile):
    f = bz2.BZ2File(knnfile, 'rb')
    knn = cPickle.load(f)
    return knn


def find(qvec, db):
    q_trans = np.array([np.insert(qvec, 0, 0)])
    distances, indices = db.kneighbors(q_trans)

    return distances, indices