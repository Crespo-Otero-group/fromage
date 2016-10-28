# Atom class to interface between formatting
# in different programs


class Atom:
        # Atom class has element type, coordinates and charge

    def __init__(self, elemIn="H", xIn=0.0, yIn=0.0, zIn=0.0, qIn=0.0):
        self.elem = elemIn
        self.x = xIn
        self.y = yIn
        self.z = zIn
        self.q = qIn

        # toString methods to be used mainly for debugging
    def __repr__(self):
        return str(self.elem) + "\t" + str(self.x) + "\t" + str(self.y) + "\t" + str(self.z) + "\t" + str(self.q)

    def __str__(self):
        return str(self.elem) + "\t" + str(self.x) + "\t" + str(self.y) + "\t" + str(self.z) + "\t" + str(self.q)
