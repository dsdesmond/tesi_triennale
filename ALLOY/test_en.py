import numpy as np
from Lattice import Lattice

a = np.ones((4, 4)); a[1, 1] = -1                 # una impurità isolata
L = Lattice(4, 1.0, "Kawasaki", a)
print(L.compute_E([1, 1]), L.compute_E([1, 2]))   # attesi: -4.0  2.0
print(-2*L.compute_E([1, 1]) - 2*L.compute_E([1, 2]) + 4*L.J)   # atteso: 0.0

b = np.ones((4, 4)); b[1, 1] = -1; b[1, 2] = -1   # due impurità vicine
M = Lattice(4, 1.0, "Kawasaki", b)
print(-2*M.compute_E([1, 2]) - 2*M.compute_E([1, 3]) + 4*M.J)   # atteso: -4.0