import numpy as np
import matplotlib.pyplot as plt

def generate_random_array(shape_of_array,list_of_peptides,number_of_blocks):
    minimal_number_of_spots = len(list_of_peptides)*number_of_blocks
    maximal_number_of_spots = shape_of_array[0]*shape_of_array[1]
    if minimal_number_of_spots > maximal_number_of_spots:
        raise ValueError('To many peptides and/blocks')




