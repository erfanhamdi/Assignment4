import os
import argparse

Nxs = [10, 40, 80]
Nys = [5, 10]
element_types = ['P', 'Q']
element_degrees = [1, 2]
quadrature_degrees = [1, 2, 4]
results_dir = "./results"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

for Nx in Nxs:
    for Ny in Nys:
        for element_type in element_types:
            for element_degree in element_degrees:
                for quadrature_degree in quadrature_degrees:
                    print(f"Running Nx={Nx}, Ny={Ny}, element_type={element_type}, element_degree={element_degree}, quadrature_degree={quadrature_degree}")
                    output_dir = f"./results/Nx_{Nx}_Ny_{Ny}_element_type_{element_type}_element_degree_{element_degree}_quadrature_degree_{quadrature_degree}"
                    # create output directory if it doesn't exist
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)
                    os.system(f"python main.py --Nx {Nx} --Ny {Ny} --element_type {element_type} --element_degree {element_degree} --quadrature_degree {quadrature_degree} --output_dir {output_dir}")




