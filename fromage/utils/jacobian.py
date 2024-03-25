import numpy as np

def update_jac(jac, mol_a, mol_b):
    for i, j in zip(mol_a, mol_b):
        jac[i, j] = 1
        jac[i + 1, j + 1] = 1
        jac[i + 2, j + 2] = 1
    return jac


def jacobian(real_atoms, aug_mol_atoms, lac_atoms, lah_atoms, la_atoms):
    """Construct M x M Jacobian from real and model regions"""

    # make real region
    M = len(real_atoms) * 3

    # get initial matrix
    jac = np.zeros((M, M))

    # get useful indices
    # inner_r = [aug_mol_atoms.index(atom) for atom in aug_mol_atoms if atom not in lac_atoms or atom not in lac_atoms]
    inner_r, lac_r, la_r = [], [], []

    for i, atom in enumerate(aug_mol_atoms):
        i *= 3
        if atom in lac_atoms:
            lac_r.append(i)
        elif atom in la_atoms:
            la_r.append(i)
        else:
            inner_r.append(i)

    inner_R, lac_R, lah_R, outer_R = [], [], [], []
    for i, atom in enumerate(real_atoms):
        i *= 3
        if atom in lac_atoms:
            lac_R.append(i)
        elif atom in lah_atoms:
            lah_R.append(i)
        elif atom in aug_mol_atoms:
            inner_R.append(i)
        else:
            outer_R.append(i)

    jac = update_jac(jac, inner_r, inner_R)
    jac = update_jac(jac, lac_r, lac_R)
    jac = update_jac(jac, outer_R, outer_R)

    for i, j in zip(la_r, lah_R):

        g = aug_mol_atoms[int(i / 3)].gfac
        jac[i, j] = g
        jac[i + 1, j + 1] = g
        jac[i + 2, j + 2] = g

    for i, j in zip(la_r, lac_R):
        # print(i)
        g = aug_mol_atoms[int(i / 3)].gfac
        jac[i, j] = 1 - g
        jac[i + 1, j + 1] = 1 - g
        jac[i + 2, j + 2] = 1 - g

    # print(inner_r, lac_r, la_r)
    # print(inner_R, lac_R, lah_R, outer_R)

    return jac

def transform_grads(in_grad,jacobian, linkatoms):
    """
    Function to transform gradients with Jacobian matrix

    First reshapes gradients to (1,N) to ensure dimensionality of grads
    matches jacobian, then mulitplies by (N,N) matrix. Returns gradients
    in same form as in_grad (N). 

    Parameters
    ----------
    in_grad : 1-D np.array
        gradients of augmented model calculation
    jac : N-D np.array
        jacobian matrix 
    
    Returns
    -------
    out_grad : 1-D np.array
        Jacobian-transformed augmented model gradients
    """


    len_LA = len(linkatoms)*3
    len_r  =len(in_grad) # length of augmented model region
    len_R = np.shape(jacobian)[0] # length of real region
    in_grad = in_grad.reshape((1,len_r)) # dim to match jacobian

    # pad augmented model gradient with zeroes to length R
    padded_grad = np.hstack([in_grad, np.zeros((1,(len_R-len_r)))])
    #print(len(padded_grad))
    # transform with Jacobian matrix
    out_grad = np.matmul(padded_grad, jacobian)
    model_len = len_r-len_LA
    out_grad = out_grad[0][:model_len]
    #Remove extra zereos
    
    print("modellen", model_len)
    in_grad = in_grad[0][:model_len].reshape(model_len)
    print("OUT_GRAD\n\n\n:", out_grad)
    print("Pre-Jacobian transform:", in_grad)
    print("Post-Jacobian transform:", out_grad)
    print("Change in gradient:", in_grad - out_grad)
    # gradients = padded_gradients - (zeroes + modified LAH grads)
    return out_grad


def jac_transform(grads, jac, linkatoms):
    """"
    projects aug_mol gradients onto model region   
    """

    #print("\n\n\nin_pos:\n", len(in_pos))
    en_gr_list = list(grads)
    init_grads = en_gr_list[1]

    en_gr_list[1] =  transform_grads(init_grads, jac,linkatoms)
    print("\n\n\nen_gr_list[1]:\n", len(en_gr_list[1]))
    en_gr_list = tuple(en_gr_list)


    return en_gr_list


if __name__ == "__main__":
    from fromage.utils.mol import Mol
    from fromage.utils.atom import Atom
    from fromage.utils.linkatom import LinkAtom

    # # set up relevent Mol objects to create Jacobian
    real_atoms = Mol(
        [
            Atom("H", -6.20823, 1.79861, -0.68813),
            Atom("O", -6.29242, 1.12771, 0.01592),
            Atom("O", -7.74237, 1.03423, -0.01592),
            Atom("H", -7.82656, 0.36332, 0.68813),
        ]
    )

    mol_atoms = Mol(
        [Atom("H", -6.20823, 1.79861, -0.68813), Atom("O", -6.29242, 1.12771, 0.01592)]
    )

    lac_atoms = Mol(Atom("O", -6.29242, 1.12771, 0.01592))
    lah_atoms = Mol(Atom("O", -7.74237, 1.03423, -0.01592))
    la_atoms = Mol(LinkAtom(lac_atoms[0], lah_atoms[0]))

    aug_mol_atoms = mol_atoms + la_atoms
    print(aug_mol_atoms, "\n", la_atoms, "\n", lac_atoms, "\n", lah_atoms)
    jac = jacobian(real_atoms, aug_mol_atoms, lac_atoms, lah_atoms, la_atoms)

    print(jac)

    grads = np.array([[0.1, 0.1, 0.1, 0.1, 0.1, 0.1]]) # h1  # o1  # h_la
    print(np.shape(grads))
    print(np.shape(jac))
    grads = transform_grads(grads, jac, la_atoms)

    print(real_atoms)