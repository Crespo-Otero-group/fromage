"""Defines the link atom object"""

from fromage.utils.atom import Atom
from fromage.utils.per_table import periodic 

class LinkAtom(Atom):
    """
    Object representing a link atom.

    Required for constructing augmented model region when covalent bonds are cut during ONIOM calculation.
    Inherits attributes and methods of Atom object. 

    Attributes
    ----------
    lac_partner : Atom
        atom in model region
    lah_partner :
        atom in real region that is being replaced by link atom
    g_fac : float
        scaling factor for position along LAC-LAH bond vector
    """
    
    def __init__(self, lac_partner, lah_partner, gfac=None, elemIn="H", xIn=0.0, yIn=0.0, zIn=0.0, qIn=0.0, num=1):
        super().__init__(elemIn, xIn,yIn,zIn, qIn,num)
        self.lac_partner = lac_partner
        self.lah_partner = lah_partner
        self.gfac = gfac
        
        # scaling factor 
        self.lac_rad = periodic[self.lac_partner.elem.lower()]["cov"]
        self.lah_rad = periodic[self.lah_partner.elem.lower()]["cov"]
        self.rad = periodic[self.elem.lower()]["cov"]
        
        if self.gfac == None:
            self.calc_gfac(self.lac_rad, self.lah_rad, self.rad)


        #initialise position
        self.initial_pos()
 
    def update_pos(self, lac_pos, lah_pos):
        """
        update LA position using new LAC/LAH coordinates

        Eqn. 19 from Chung 2015 (https://pubs.acs.org/doi/pdf/10.1021/cr5004419)
        """
        self.set_pos(lac_pos+self.gfac*(lah_pos-lac_pos))
        return 
    
    def initial_pos(self):
        """for first instance of class"""
        self.update_pos(self.lac_partner.get_pos(), self.lah_partner.get_pos())
        return 
    
    def calc_gfac(self, lac_rad, lah_rad, rad):
        """Calculate scaling factor using covalent radii"""
        self.gfac = (lac_rad + rad)/(lah_rad + lac_rad) 
        return 

    def set_gfac(self, new_gfac):
        """set custom scaling factor"""
        self.gfac = new_gfac
        return

if __name__=='__main__':
    #illustrating functionality 
    lac_atom = Atom("C", 0,0,0)
    lah_atom = Atom("C", 1, 1, 1)
    
    la_atom = LinkAtom(lac_atom,lah_atom, gfac=0.75)
    
    

    print("LAC atom: ", lac_atom)

    print("LAH atom: ", lah_atom)
    print("New link atom atom: ", la_atom)
    
