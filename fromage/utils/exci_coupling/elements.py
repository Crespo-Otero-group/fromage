""" Shortened periodic table data"""

class element():
    """
    Class representing an element with information on the symbol,name,
    atomic number and relative atomic mass

    Attributes
    ----------
    symbol : String
        Atomic symbol is passed to the class
    name: String
        Element name
    mass: float
        Relative atomic mass associated with the passed symbol
    number: integer
        Atomic number associated with the passed symbol
    """
    def __init__(self,identifier):

       names={"H": "Hydrogen","B":"Boron","C":"Carbon","N":"Nitrogen","O":"Oxygen",
       "F":"Fluorine","Si":"Silicon","P":"Phosphorus","S":"Sulfur","Cl":"Chlorine",
             "Br":"Bromine","I":"Iodine"}
       masses={"H": 1.007825, "B":10.81,"C": 12.011, "N":14.007, "O":15.999,
              "F":18.998,"Si":28.0855,"P":30.9738,"S":32.06,"Cl":35.453,"Br":79.904,
              "I":126.909}
       symb_to_numbers={"H": 1, "B":5, "C": 6, "N": 7, "O": 8, "F":9,"Si":14,"P":15,"S":16,
               "Cl":17,"Br":35,"I":53}
       numbers_to_symb={ 1:"H", 5:"B", 6:"C", 7:"N", 8:"O", 9:"F",14:"Si",15:"P",16:"S",
               17:"Cl",35:"Br",53:"I"}

       if isinstance(identifier, int):
           self.symbol = numbers_to_symb[identifier]
           self.name = names[self.symbol]
           self.mass = masses[self.symbol]
           self.atomic = identifier

       elif isinstance(identifier, str):
           self.symbol = identifier
           self.name = names[identifier]
           self.mass = masses[identifier]
           self.atomic = symb_to_numbers[identifier]
