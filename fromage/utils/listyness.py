"""Useful list functions to import

This is handy to make lists (or lists of lists of-, and so on) of Atoms into
classes. We want list-like behaviour, be we do not want to inherit max() or
min(), or any other number of senseless functions. This is also cleaner than
an abstract class, because we can easily change the chief_list attribute,
making it clearer for the user.

"""

def append(self, element):
    getattr(self,self.chief_list).append(element)


def extend(self, other):
    getattr(self,self.chief_list).extend(getattr(other,other.chief_list))


def insert(self, i, element):
    getattr(self,self.chief_list).insert(i, element)


def remove(self, element):
    getattr(self,self.chief_list).remove(element)


def index(self, element):
    return getattr(self,self.chief_list).index(element)


def pop(self, i=-1):
    return getattr(self,self.chief_list).Pop(i)


def clear(self):
    getattr(self,self.chief_list).clear()


def count(self, element):
    return getattr(self,self.chief_list).count()

def __len__(self):
    return len(getattr(self,self.chief_list))


def __eq__(self, other):
    return getattr(self,self.chief_list) == other.self.chief_list


def __getitem__(self, index):
    return getattr(self,self.chief_list)[index]


def __setitem__(self, index, value):
    getattr(self,self.chief_list)[index] = value
    return


def __contains__(self, elem):
    return getattr(self,self.chief_list).__contains__(elem)
