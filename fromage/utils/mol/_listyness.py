from copy import deepcopy

# list-like behaviour


def append(self, element):
    self.atoms.append(element)


def extend(self, other_mol):
    self.atoms.extend(other_mol.atoms)


def insert(self, i, element):
    self.atoms.insert(i, element)


def remove(self, element):
    self.atoms.remove(element)


def index(self, element):
    return self.atoms.index(element)


def pop(self, i=-1):
    return self.atoms.pop(i)


def clear(self):
    self.atoms.clear()


def count(self, element):
    return self.atoms.count()


def __add__(self, other_mol):
    import fromage.utils.mol as mol_init
    mol_init.try_ismol(other_mol)
    return mol_init.Mol(deepcopy(self).atoms + other_mol.atoms, vectors=self.vectors, bonding=self.bonding, thresh=self.thresh)


def __len__(self):
    return len(self.atoms)


def __eq__(self, other):
    return self.atoms == other.atoms


def __getitem__(self, index):
    return self.atoms[index]


def __setitem__(self, index, value):
    self.atoms[index] = value
    return


def __contains__(self, elem):
    return self.atoms.__contains__(elem)
