import numpy as np
from ufl import sym, grad, Identity, tr, inner, Measure, TestFunction, TrialFunction

from mpi4py import MPI

from dolfinx import fem, io
import dolfinx.fem.petsc
from dolfinx.mesh import create_rectangle, CellType

import argparse
import os
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("--Nx", type=int, default=80)
arg_parser.add_argument("--Ny", type=int, default=10)
arg_parser.add_argument("--element_type", type=str, default='P')
arg_parser.add_argument("--element_degree", type=int, default=2)
arg_parser.add_argument("--quadrature_degree", type=int, default=2)
arg_parser.add_argument("--output_dir", type=str, default="/projectnb/me700/students/erfan/4-2/results")

args = arg_parser.parse_args()
if args.element_type == 'P':
    cell_type = CellType.triangle
elif args.element_type == 'Q':
    cell_type = CellType.quadrilateral

length, height = 10, 1.0
Nx, Ny = args.Nx, args.Ny
output_dir = args.output_dir
domain = create_rectangle(
    MPI.COMM_WORLD,
    [np.array([0, 0]), np.array([length, height])],
    [Nx, Ny],
    cell_type=cell_type,
)

dim = domain.topology.dim

degree = args.element_degree
shape = (dim,)  # this means we want a vector field of size `dim`
V = fem.functionspace(domain, ("Lagrange", degree, shape))

num_dofs_global = V.dofmap.index_map.size_global * V.dofmap.index_map_bs

u_sol = fem.Function(V, name="Displacement")

E = fem.Constant(domain, 210e3)
nu = fem.Constant(domain, 0.3)
I = fem.Constant(domain, (1/12) * 1 * height**3)
lmbda = E * nu / (1 + nu) / (1 - 2 * nu)
mu = E / 2 / (1 + nu)


def epsilon(v):
    return sym(grad(v))


def sigma(v):
    return lmbda * tr(epsilon(v)) * Identity(dim) + 2 * mu * epsilon(v)

u = TrialFunction(V)
v = TestFunction(V)

rho = 2e-3
g = 9.81
f = fem.Constant(domain, np.array([0, -rho * g]))

dx = Measure("dx", domain=domain, metadata={"quadrature_degree": args.quadrature_degree})
a = inner(sigma(u), epsilon(v)) * dx
L = inner(f, v) * dx

def left(x):
    return np.isclose(x[0], 0)


def right(x):
    return np.isclose(x[0], length)


left_dofs = fem.locate_dofs_geometrical(V, left)
right_dofs = fem.locate_dofs_geometrical(V, right)

V_uy, mapping = V.sub(1).collapse()
right_dofs_uy = fem.locate_dofs_geometrical((V.sub(1), V_uy), right)
left_dofs_uy = fem.locate_dofs_geometrical((V.sub(1), V_uy), left)
uD_y = fem.Function(V_uy)
bcs = [
    fem.dirichletbc(uD_y, left_dofs_uy, V.sub(1)),
    fem.dirichletbc(uD_y, right_dofs_uy, V.sub(1)),
]

problem = fem.petsc.LinearProblem(
    a, L, u=u_sol, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"}
)
def monitor(ksp, its, rnorm):
    print(f"Iteration {its}, residual {rnorm}")

problem.solver.setMonitor(monitor)

problem.solve()
if not problem.solver.is_converged:
    print("Solver did not converge")
    exit(1)

tol = 0.001  # Avoid hitting the outside of the domain
x = np.linspace(0 + tol, length - tol, 101)
points = np.zeros((3, 101))
points[0] = x
points[1] = height/2
u_values = []
p_values = []

from dolfinx import geometry
bb_tree = geometry.bb_tree(domain, domain.topology.dim)

cells = []
points_on_proc = []
# Find cells whose bounding-box collide with the the points
cell_candidates = geometry.compute_collisions_points(bb_tree, points.T)
# Choose one of the cells that contains the point
colliding_cells = geometry.compute_colliding_cells(domain, cell_candidates, points.T)
for i, point in enumerate(points.T):
    if len(colliding_cells.links(i)) > 0:
        points_on_proc.append(point)
        cells.append(colliding_cells.links(i)[0])

points_on_proc = np.array(points_on_proc, dtype=np.float64)
u_values = u_sol.eval(points_on_proc, cells)
np.save(os.path.join(output_dir, f"dofs_{num_dofs_global}_u_values.npy"), u_values)

import matplotlib.pyplot as plt
fig = plt.figure()
plt.plot(points_on_proc[:, 0], u_values[:, 1], "k", linewidth=2, label="FEA")
plt.grid(True)

u_analytic = lambda x: (-rho*g/(24*E.value*I.value))*(x**4 - 2*length*x**3 + length**3*x)
x_analytic = np.linspace(0 + tol, length - tol, 101)
u_analytic_values = u_analytic(x_analytic)

error = np.max(np.abs(u_values[:, 1] - u_analytic_values))
np.savetxt(os.path.join(output_dir, f"error.txt"), [error])
print(f"Error: {error}")
plt.plot(x_analytic, u_analytic_values, "r", linewidth=2, label="Analytic")
plt.grid(True)
plt.xlabel("x")
plt.ylabel("displacement")
plt.legend()
output_png = os.path.join(output_dir, f"deflection.png")
output_pvd = os.path.join(output_dir, f"linear_elasticity.pvd")
plt.savefig(output_png)
# vtk = io.VTKFile(domain.comm, output_pvd, "w")
# vtk.write_function(u_sol)
# vtk.close()