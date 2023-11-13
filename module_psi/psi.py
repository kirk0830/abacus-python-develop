import numpy as np

class psi:
    
    current_k = 0
    current_b = 0

    def __init__(self, nkpts: int, nbands: int, npwx: int):
        self.nkpts = nkpts
        self.nbands = nbands
        self.npwx = npwx
        self.psi = np.zeros((nkpts, nbands, npwx), dtype=np.complex128)

    def fix_k(self, k: int):
        self.current_k = k

    def fix_b(self, b: int):
        self.current_b = b

    def get_nbands(self):
        return self.nbands
    
    def get_nkpts(self):
        return self.nkpts
    
    def get_nbasis(self):
        return self.npwx
    
    def get_psi(self):
        return self.psi
    
    def v(self, *args):
        if len(args) == 1:
            return self.psi[self.current_k, self.current_b, args[0]]
        elif len(args) == 2:
            return self.psi[self.current_k, args[0], args[1]]
        elif len(args) == 3:
            return self.psi[args[0], args[1], args[2]]
        else:
            raise ValueError("Wrong number of arguments")
        
    def print(self):
        for ikpt in range(self.nkpts):
            for iband in range(self.nbands):
                print(f'psi[{ikpt+1}, {iband+1}, :] = {self.psi[ikpt, iband, :]}')

    def print_k(self, k: int):
        for iband in range(self.nbands):
            print(f'psi[{k+1}, {iband+1}, :] = {self.psi[k, iband, :]}')

    def print_b(self, k: int, b: int):
        print(f'psi[{k+1}, {b+1}, :] = {self.psi[k, b, :]}')

    def overlap_with(self, ikpt, other):

        return self.sandwich_with(
            ikpt, 
            np.ones((self.npwx, self.npwx), dtype=np.complex128), 
            other)

    def sandwich_with(self, ikpt, operator, other):
        if self.nkpts != other.nkpts:
            raise ValueError("Number of k-points does not match")
        elif self.nbands != other.nbands:
            raise ValueError("Number of bands does not match")
        elif self.npwx != other.npwx:
            raise ValueError("Number of orbitals does not match")
        else:
            sandwich = np.zeros((self.nbands, self.nbands), dtype=np.complex128)
            for iband in range(self.nbands):
                for jband in range(self.nbands):
                    sandwich[iband, jband] = np.vdot(self.psi[ikpt, iband, :], operator @ other.psi[ikpt, jband, :])
            return sandwich

    def calculate_norm(self):
        norm = np.zeros((self.nkpts, self.nbands), dtype=np.float64)
        for ikpt in range(self.nkpts):
            for iband in range(self.nbands):
                norm[ikpt, iband] = np.linalg.norm(self.psi[ikpt, iband, :])
        return norm
    