import numpy as np
a = np.array([[[1, 2], [3, 4]],[[5, 6], [7, 8]]])

b = np.array([[[1, 2], [3, 4]],[[5, 6], [7, 8]]])
d = np.array([[[1, 2], [3, 4]],[[5, 6], [7, 8]]])
c = np.concatenate((a, b, d), axis=2)
print(c)
print(c.shape)


v_min = c.min(axis=(0, 1), keepdims=True)
v_max = c.max(axis=(0, 1), keepdims=True)

print(v_min)
print(v_max)
norm = (c - v_min)/(v_max - v_min)
print(norm)
