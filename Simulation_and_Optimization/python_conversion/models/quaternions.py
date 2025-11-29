
# small list of quaternion functions

import numpy as np

# this is hamilton's quaternion product
def product(a, b):
    v = b[0] * a[1:] + a[0] * b[1:] + np.cross(a[1:], b[1:])
    return np.array([a[0] * b[0] - np.dot(a[1:], b[1:]), v[0], v[1], v[2]])

# this is the inverse of a unit quat
def conjugate(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])

# this rotates a vector with a fixed frame
def sandwich(q, v):
    return product(q, product(np.array([0, v[0], v[1], v[2]]), conjugate(q)))[1:]

# this rotates a frame with a fixed vector
def frame_rotation(q, v):
    return sandwich(conjugate(q), v)

# this constrains a quat to S3 or vector to S2
def normalize(q):
    norm = np.linalg.norm(q)
    norm = norm if norm !=0 else 1
    return q / norm

