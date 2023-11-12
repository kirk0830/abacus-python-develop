import numpy as np

def generate_pw_basis(ecutwfc: float, unitcell: dict):

    """
    Generate a plane wave basis set for a given cutoff energy and cell dimensions.
    """
    g_vectors = []
    nbx = int(np.sqrt(ecutwfc * unitcell['GGT'][0, 0]) / unitcell['tpiba']) + 1
    nby = int(np.sqrt(ecutwfc * unitcell['GGT'][1, 1]) / unitcell['tpiba']) + 1
    nbz = int(np.sqrt(ecutwfc * unitcell['GGT'][2, 2]) / unitcell['tpiba']) + 1
    for ibx in range(-nbx, nbx + 1):
        for iby in range(-nby, nby + 1):
            for ibz in range(-nbz, nbz + 1):
                if ibx**2 * unitcell['GGT'][0, 0] + iby**2 * unitcell['GGT'][1, 1] + ibz**2 * unitcell['GGT'][2, 2] < ecutwfc:
                    # collect the G-vectors
                    g_vectors.append([ibx, iby, ibz])
    g_vectors = np.array(g_vectors)
    print("PW Basis generation done, number of G-vectors: ", len(g_vectors))
    return g_vectors
    