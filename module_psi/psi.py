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