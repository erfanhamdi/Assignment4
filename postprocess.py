import numpy as np
import matplotlib.pyplot as plt
import glob
def u_analytic():
    length = 10
    height = 1
    E = 210e3
    nu = 0.3
    I = (1/12) * 1 * height**3
    rho  = 2e-3
    g = 9.81
    tol = 0.001
    x = np.linspace(0 + tol, length - tol, 101)
    u_analytic = (-rho*g/(24*E*I))*(x**4 - 2*length*x**3 + length**3*x)
    return u_analytic
def main():
    plt.figure()
    u_analytic_ = u_analytic()
    element_types = ['P', 'Q']
    element_degrees = [1, 2]
    quadrature_degrees = [2]
    error_dict = {}
    for element_type in element_types:
        for element_degree in element_degrees:
            for quadrature_degree in quadrature_degrees:
                results = glob.glob(f"./results/*element_type_{element_type}_element_degree_{element_degree}_quadrature_degree_{quadrature_degree}/dofs_*u_values.npy")
                # sort the results by dofs
                results.sort(key=lambda x: int(x.split("/")[-1].split("_")[1]))
                errors_list = []
                dofs_list = []
                for result in results:
                    u_values = np.load(result)
                    dofs = result.split("/")[-1].split("_")[1]
                    dofs_list.append(int(dofs))
                    element_type = result.split("/")[-2].split("_")[6]
                    element_degree = result.split("/")[-2].split("_")[9]
                    quadrature_degree = result.split("/")[-2].split("_")[12]
                    print(f"dofs: {dofs}, element_type: {element_type}, element_degree: {element_degree}, quadrature_degree: {quadrature_degree}")
                    error = np.max(np.abs(u_values[:, 1] - u_analytic_))
                    errors_list.append(error)
                    out_array = np.array([dofs_list, errors_list])
                error_dict[f"{element_type}_{element_degree}_{quadrature_degree}"] = out_array
    line_styles = ['-', '--', '-.', ':']
    for i, element in enumerate(error_dict.keys()):
        dofs = list(error_dict[element][0, :])
        errors = list(error_dict[element][1, :])
        plt.plot(dofs, errors, label=f"{element}", linestyle=line_styles[i])
        # plt.scatter(dofs, errors, label=f"{element}")
        plt.xlabel("#dofs")
        plt.ylabel("Max Error")
    plt.legend()
    plt.savefig("errors.png")
    plt.show()
if __name__ == "__main__":
    main()
