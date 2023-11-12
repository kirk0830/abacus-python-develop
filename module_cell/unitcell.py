import numpy as np

def read_atomic_species(file_handle, line: str, clean_line: str) -> dict:
    atomic_species = {}
    
    while clean_line != "\n":
        element_symbol = clean_line.split()[0]
        mass = float(clean_line.split()[1])
        pseudo_file = clean_line.split()[2]
        atomic_species[element_symbol] = {"mass": mass, "pseudo_file": pseudo_file}
        line = file_handle.readline()
        clean_line = line.strip()

    return atomic_species

def read_numerical_orbitals(file_handle, line: str, clean_line: str) -> dict:
    numerical_orbitals = {}
    while clean_line != "\n":
        words = clean_line.split()
        numerical_orbitals[words[0]] = clean_line
        line = file_handle.readline()
        clean_line = line.strip()
    return numerical_orbitals

def read_lattice_constant(file_handle, line: str, clean_line: str) -> float:
    lattice_constant = float(clean_line.split()[0])
    return lattice_constant

def read_lattice_vectors(file_handle, line: str, clean_line: str) -> np.ndarray:
    lattice_vectors = np.zeros((3, 3))
    for i in range(3):
        if i > 0:
            line = file_handle.readline()
            clean_line = line.strip()
        lattice_vectors[i, 0] = float(clean_line.split()[0])
        lattice_vectors[i, 1] = float(clean_line.split()[1])
        lattice_vectors[i, 2] = float(clean_line.split()[2])
    return lattice_vectors

def read_atomic_positions(file_handle, line: str, clean_line: str, ntypes: int) -> dict:
    atomic_positions = {}
    number_of_atomtypes_read = 0
    while number_of_atomtypes_read < ntypes:
        element_symbol = clean_line.split()[0]
        atomic_positions[element_symbol] = []
        line = file_handle.readline()
        clean_line = line.strip()
        number_of_atoms = int(clean_line.split()[0])
        atomic_positions[element_symbol].append(number_of_atoms)
        line = file_handle.readline()
        clean_line = line.strip()
        starting_magnetization = float(clean_line.split()[0])
        atomic_positions[element_symbol].append(starting_magnetization)
        number_of_atoms_read = 0
        while number_of_atoms_read < number_of_atoms:
            line = file_handle.readline()
            clean_line = line.strip()
            x = float(clean_line.split()[0])
            y = float(clean_line.split()[1])
            z = float(clean_line.split()[2])
            atomic_positions[element_symbol].append([x, y, z])
            number_of_atoms_read += 1
        number_of_atomtypes_read += 1
    return atomic_positions

def import_structure(file: str) -> dict:
    structure = {}

    b_read_atomic_species = False
    b_read_numerical_orbitals = False
    b_read_lattice_constant = False
    b_read_lattice_vectors = False
    b_read_atomic_positions = False
    line = "start"
    with open(file, 'r') as f:
        while line:
            line = f.readline()
            clean_line = line.strip()

            if clean_line.startswith("ATOMIC_SPECIES"):
                b_read_atomic_species = True
            if clean_line.startswith("NUMERICAL_ORBITALS"):
                b_read_numerical_orbitals = True
            if clean_line.startswith("CELL_CONSTANT"):
                b_read_lattice_constant = True
            if clean_line.startswith("CELL_PARAMETERS"):
                b_read_lattice_vectors = True
            if clean_line.startswith("ATOMIC_POSITIONS"):
                b_read_atomic_positions = True

            if b_read_atomic_species:
                structure["atomic_species"] = read_atomic_species(f, line, clean_line)
                b_read_atomic_species = False
            if b_read_numerical_orbitals:
                structure["numerical_orbitals"] = read_numerical_orbitals(f, line, clean_line)
                b_read_numerical_orbitals = False
            if b_read_lattice_constant:
                structure["lattice_constant"] = read_lattice_constant(f, line, clean_line)
                b_read_lattice_constant = False
            if b_read_lattice_vectors:
                structure["lattice_vectors"] = read_lattice_vectors(f, line, clean_line)
                b_read_lattice_vectors = False
            if b_read_atomic_positions:
                structure["atomic_positions"] = read_atomic_positions(f, line, clean_line, len(structure["atomic_species"]))
                b_read_atomic_positions = False

    print("Structure import done. Number of atom types: ", len(structure["atomic_species"]))
    return structure

def generate_unitcell(
        lattice_constant: float, 
        a: float, b: float, c:float, 
        alpha: float, beta: float, gamma: float,
        structure: dict) -> dict:
    """
    Generate a unit cell from the lattice parameters a, b, c, alpha, beta, gamma.
    """
    cell = np.zeros((3, 3))
    cell[0, 0] = a
    cell[1, 0] = b * np.cos(gamma)
    cell[1, 1] = b * np.sin(gamma)
    cell[2, 0] = c * np.cos(beta)
    cell[2, 1] = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cell[2, 2] = np.sqrt(c**2 - cell[2, 0]**2 - cell[2, 1]**2)
    
    tpiba = 2 * np.pi / lattice_constant

    GT = np.linalg.inv(cell)
    G = GT.T
    GGT = np.dot(G, GT)
    
    unitcell = {}
    unitcell["cell"] = {
        "a": a,
        "b": b,
        "c": c,
        "alpha": alpha,
        "beta": beta,
        "gamma": gamma,
        "latvec": cell,
        "GT": GT,
        "G": G,
        "GGT": GGT,
        "tpiba": tpiba
    }
    unitcell["structure"] = structure

    return unitcell
