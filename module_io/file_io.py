def import_numerical_atomic_orbitals(file: str) -> dict:
    numerical_atomic_orbitals = {}

    line = "start"
    with open(file, 'r') as f:
        while line:
            line = f.readline()
            clean_line = line.strip()
            if clean_line.startswith("Type"):
                # then it is a new numerical orbital
                line = f.readline()
                clean_line = line.strip()
                type_id = clean_line.split()[0]
                l = clean_line.split()[1]
                n = clean_line.split()[2]

    return numerical_atomic_orbitals

def import_lowf_file(file: str) -> dict:
    band_structure = {}

    line = "start"
    with open(file, 'r') as f:
        while line:
            pass

    return band_structure