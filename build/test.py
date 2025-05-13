import numpy as np
import matplotlib.pyplot as plt
from adios2 import Adios, Stream, bindings


num_points = 41
filename = "data.bp"
domain_size = 1.0

adios = Adios()
io = adios.declare_io("SimulationOutput")
stream = Stream(io, filename, "r")

u = v = p = None


while True:
    status = stream.begin_step()
    if status != bindings.StepStatus.OK:
        break

    var_u = stream.inquire_variable("u")
    var_v = stream.inquire_variable("v")
    var_p = stream.inquire_variable("p")

    if var_u and var_v and var_p:
        u = stream.read(var_u)
        v = stream.read(var_v)
        p = stream.read(var_p)
    else:
        print("Error: Could not find all required variables.")
        break

    stream.end_step()

stream.close()

x = np.linspace(0, domain_size, num_points)
y = np.linspace(0, domain_size, num_points)
X, Y = np.meshgrid(x, y)

u_next = np.reshape(u, (num_points, num_points))
v_next = np.reshape(v, (num_points, num_points))
p_next = np.reshape(p, (num_points, num_points))

plt.style.use("dark_background")
plt.figure(figsize=(8, 6))


contour = plt.contourf(X[::2, ::2], Y[::2, ::2], p_next[::2, ::2], cmap="coolwarm")
plt.colorbar(contour)

plt.quiver(X[::2, ::2], Y[::2, ::2], u_next[::2, ::2], v_next[::2, ::2], color="black")


plt.xlim((0, 1))
plt.ylim((0, 1))
plt.xlabel("x")
plt.ylabel("y")
plt.title("Final Timestep: Pressure Contour and Velocity Field")

plt.tight_layout()
plt.show()
