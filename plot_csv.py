import numpy as np
import matplotlib.pyplot as plt

with open("sundials2_sol.csv") as file_name:
    t, x, y = np.loadtxt(file_name, delimiter=",", comments="#", unpack=True)

plt.plot(x,y)
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("xysolution.png", dpi=400, bbox_inches="tight")
plt.show()

plt.plot(t, x, label="x")
plt.plot(t, y, label="y")
plt.xlabel("t")
plt.ylabel("solutions")
plt.legend()
plt.show()

