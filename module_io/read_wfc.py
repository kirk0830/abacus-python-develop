import module_psi.psi as psi
import numpy as np

def read_lowf(path: str, nkpts: int, nbands_ref: int, nbasis_ref: int) -> psi.psi:

    if path.endswith('/'):
        path = path[:-1]

    for ikpt in range(nkpts):
        nbands = 0
        nbasis = 0
        band_energies = []
        band_occupations = []
        coefficients = ''
        with open(path + f'/LOWF_K_{ikpt+1}.txt', 'r') as f:
            lines = f.readlines()
        for line in lines:
            # remove /n, /t, /r and spaces at the beginning and end of the line
            line = line.strip()
            # skip empty lines
            if len(line) == 0:
                continue
            
            if line.endswith('(index of k points)'):
                continue
            elif line.endswith('(number of bands)'):
                nbands = int(line.split()[0])
                continue
            elif line.endswith('(number of orbitals)'):
                nbasis = int(line.split()[0])
                continue
            elif line.endswith('(band)'):
                continue
            elif line.endswith('(Ry)'):
                continue
            elif line.endswith('(Occupations)'):
                continue
            else:
                if len(line.split())%2 != 0:
                    continue
                else:
                    coefficients += line + ' '
                    continue
        coefficients = coefficients.split()


        if len(coefficients) != nbands * nbasis * 2:
            raise ValueError(f'For kpoint ({ikpt}), number of coefficients ({len(coefficients)}) does not match the number of bands ({nbands}) and orbitals ({nbasis})')
        elif nbands != nbands_ref or nbasis != nbasis_ref:
            raise ValueError(f'Number of bands ({nbands}) or orbitals ({nbasis}) does not match the reference values ({nbands_ref}, {nbasis_ref})')
        else:
            coefficients = [complex(float(coefficients[i]), float(coefficients[i+1])) for i in range(0, len(coefficients), 2)]
            coefficients = np.array(coefficients).reshape((nbands, nbasis))
            band_energies = np.array(band_energies)
            band_occupations = np.array(band_occupations)
            if ikpt == 0:
                psi_out = psi.psi(nkpts, nbands, nbasis)
            psi_out.psi[ikpt, :, :] = coefficients

    return psi_out