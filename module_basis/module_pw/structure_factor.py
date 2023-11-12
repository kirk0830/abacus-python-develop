import numpy as np

def calculate_structure_factor(g_vector: np.array, type_id: int, structure: dict) -> float:

    """
    Calculate the structure factor for a given g-vector and atom type.
    """
    structure_factor = 0
    for atom in structure["atomic_positions"][type_id]:
        structure_factor += np.exp(-2j * np.pi * np.dot(g_vector, atom))
    return structure_factor
