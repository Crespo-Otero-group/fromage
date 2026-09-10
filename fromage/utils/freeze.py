"""
Testing the freeing coordinates (just for cartesian)


for frozen atoms and frozen dihedrals, passed to the scipy optimise minimise function

> frozen atoms we just zero the gradient componenets
> frozen dihedrals, the direction along which the graident would change an angle is projected out and holds the dihedral at its starting value hopefully

"""

import numpy as np

def _as_int_list(raw):
    """ read_config to interpret frozen_at 12 to freeze say atom 12 and not 1 and 2
    etc.
    """
    if raw is None:
        return []
    if isinstance(raw, str):
        raw = raw.split()
    return [int(f) for f in raw]

def parse_frozen_atoms(raw):
    return sorted({i - 1 for i in _as_int_list(raw)})

def parse_frozen_dihedrals(raw):
    idx = _as_int_list(raw)
    if len(idx) % 4 != 0:
        raise ValueError("must have atom indices in multiples of four")
    return [tuple(i - 1 for i in idx[n:n +4]) for n in range(0, len(idx), 4)]

def dihedral_gradient(pos, quad):
    """
    pos is flattened cartesian coordinates of the otpimised region
    quad is tuple of four indices for the angle
    """

    i, j, k, l  = quad
    r = np.asarray(pos, dtype=float).reshape(-1, 3)

    f = r[i] - r[j]
    g = r[j] - r[k]
    h = r[l] - r[k]

    a = np.cross(f, g)
    b = np.cross(h, g)
    a2 = np.dot(a, a)
    b2 = np.dot(b, b)
    g_norm = np.linalg.norm(g)

    grad = np.zeros_like(r)

    if a2 < 1e-12 or b2 < 1e-12 or g_norm < 1e-12:
        return grad.ravel()


    g2 = g_norm * g_norm
    t_i = -(g_norm / a2) * a
    t_l = (g_norm / b2) * b
    c_f = np.dot(f, g) / g2
    c_h = np.dot(h, g) / g2


    #calculate the gradients
    grad[i] = -t_i
    grad[l] = -t_l
    grad[j] = (1 +c_f) * t_i + c_h * t_l
    grad[k] = -c_f * t_i + (1 - c_h) * t_l

    return grad.ravel()



def project_out_dihedral(pos, grad, quad):
    """
    Removes gradient componenet that would change one frozen idhedral

    This projects the gradient onto the subspace orhtogonal to the idhedrals 
    cartesian derivative
    """
    grad = np.array(grad, dtype=float)
    vec = dihedral_gradient(pos, quad)
    norm = np.linalg.norm(vec)
    if norm < 1e-10:
        return grad
    u = vec/norm
    return grad - np.dot(grad, u) * u

def apply_constraints(pos, grad, frozen_atoms=None, frozen_dihedrals=None):
    """
    Apply the dihedral proejction and atom freezing to a cartesian gradient
    """

    grad_arr = np.asarray(grad, dtype=float)
    if grad_arr.ndim == 0:
        return grad
    if frozen_dihedrals:
        grad_arr = project_out_dihedral(pos, grad_arr, frozen_dihedrals[0])

    if frozen_atoms:
        grad_arr = np.array(grad_arr, dtype=float)
        for i in frozen_atoms:
            grad_arr[3 * i:3 * i + 3] = 0.0

    return grad_arr
