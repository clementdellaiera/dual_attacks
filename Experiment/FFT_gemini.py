import cmath # Pour les nombres complexes (exp, etc.)
from math import pi

def addition_vecteurs(v1, v2, n):
    """
    Addition de deux vecteurs (tuples) modulo n.
    """
    return tuple((x + y) % n for x, y in zip(v1, v2))

def soustraction_vecteurs(v1, v2, n):
    """
    Soustraction de deux vecteurs (tuples) modulo n.
    """
    return tuple((x - y) % n for x, y in zip(v1, v2))

def scalar_multiplication(v, s, n):
    """
    Multiplication d'un vecteur par un scalaire s modulo n.
    """
    return tuple((x * s) % n for x in v)

def dot_product(v1, v2, n):
    """
    Produit scalaire de deux vecteurs (tuples) modulo n.
    """
    return sum(x * y for x, y in zip(v1, v2)) % n

def fft_1d_sparse(sparse_array_1d, n):
    """
    Calcule la FFT 1D d'un tableau creux.
    Entrée : un dictionnaire {indice: valeur_complexe}.
    Sortie : un nouveau dictionnaire pour la transformée.
    """
    output_fft = {}
    
    # Le nouveau tableau aura n éléments,
    # mais on ne les calcule que s'ils sont potentiellement non nuls.
    
    # Itération sur les fréquences de sortie (k de 0 à n-1)
    for k in range(n):
        sum_val = 0
        
        # Sommation sur tous les indices d'entrée i
        for i, val in sparse_array_1d.items():
            # Formule de la DFT : X_k = sum_{i=0}^{n-1} x_i * exp(-2*pi*j*i*k / n)
            angle = -2 * pi * i * k / n
            sum_val += val * cmath.exp(angle * 1j)
            
        # On ne stocke la valeur que si elle est non nulle (tolérance)
        if abs(sum_val) > 1e-9:
            output_fft[k] = sum_val
            
    return output_fft

def fft_1d_sparse_efficient(sparse_array_1d, n):
    """
    Calcule la FFT 1D d'un tableau creux en utilisant un algorithme récursif.
    """
    # Si le tableau a un seul élément
    if n == 1:
        return {0: sparse_array_1d.get(0, 0)}

    # Si le nombre de points n'est pas une puissance de 2
    if n & (n - 1) != 0:
        raise ValueError("n doit être une puissance de 2")

    # Diviser le tableau en parties paires et impaires
    even_part = {}
    odd_part = {}
    for i, val in sparse_array_1d.items():
        if i % 2 == 0:
            even_part[i // 2] = val
        else:
            odd_part[(i - 1) // 2] = val

    # Appels récursifs sur les sous-problèmes
    even_fft = fft_1d_sparse_efficient(even_part, n // 2)
    odd_fft = fft_1d_sparse_efficient(odd_part, n // 2)

    output = {}
    for k in range(n // 2):
        omega = cmath.exp(-2 * pi * k * 1j / n)
        
        # Utilisation de .get pour gérer les clés manquantes (valeur 0)
        output[k] = even_fft.get(k, 0) + omega * odd_fft.get(k, 0)
        output[k + n // 2] = even_fft.get(k, 0) - omega * odd_fft.get(k, 0)
    
    # On pourrait ajouter une étape de nettoyage pour enlever les valeurs proches de 0
    return output

def fft_multidim_sparse(sparse_table, n, k):
    """
    Calcule la FFT multi-dimensionnelle d'un tableau creux.
    
    Args:
        sparse_table (dict): Dictionnaire {(i_1, ..., i_k): valeur_complexe}.
        n (int): La taille de chaque dimension (n = 2**L).
        k (int): La dimension du vecteur d'indices.
        
    Returns:
        dict: Le dictionnaire du tableau de la transformée.
    """
    current_table = sparse_table.copy()
    
    # On fait la FFT sur chaque dimension, une par une
    for dim_index in range(k):
        next_table = {}
        
        # On regroupe les données pour faire la FFT sur la dimension courante
        # Cette étape est cruciale pour l'efficacité
        sub_tables = {}
        for indices, value in current_table.items():
            # On sépare l'indice de la dimension courante des autres
            current_dim_idx = indices[dim_index]
            other_indices_tuple = indices[:dim_index] + indices[dim_index+1:]
            
            if other_indices_tuple not in sub_tables:
                sub_tables[other_indices_tuple] = {}
            
            sub_tables[other_indices_tuple][current_dim_idx] = value
        
        # On applique la FFT 1D sur chaque sous-tableau
        for other_indices, sub_table_1d in sub_tables.items():
            # Calcule la FFT 1D creuse sur le sous-tableau
            transformed_sub_table_1d = fft_1d_sparse_efficient(sub_table_1d, n)
            
            # Reconstruit le tableau multi-dimensionnel avec les nouvelles fréquences
            for freq_1d, transformed_value in transformed_sub_table_1d.items():
                new_indices = other_indices[:dim_index] + (freq_1d,) + other_indices[dim_index:]
                if abs(transformed_value) > 1e-9: # tolérance pour la non-nullité
                    next_table[new_indices] = transformed_value
                    
        current_table = next_table
        
    return current_table

# Exemple d'utilisation
# k = 3, n = 4
k = 7
n = 8 # Doit être une puissance de 2 pour la FFT Cooley-Tukey

# Tableau d'entrée creux
T_input = {
    (1,0,0,0, 0, 0,0,0): 1.0,
    }
'''
    (1, 0, 0): 2.0,
    (0, 1, 0): 3.0,
    (0, 0, 1): 4.0,
    (2, 2, 2): 5.0,
'''

'''
# Calcule la FFT
from time import time
start = time()
T_fft = fft_multidim_sparse(T_input, n, k)
print('multidimensional FFT : ')
print('Parameters : modulo ',n,' | dimension : ', k )
print('Execution time : ', time() - start,' s')


print("Tableau de la transformée de Fourier :")
for indices, value in T_fft.items():
    print(f"{indices}: {value:.2f}")
    '''