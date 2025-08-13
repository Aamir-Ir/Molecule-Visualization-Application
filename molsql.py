#SQL.py
#import molecule
#import MolDisplay
import os
import sqlite3

class Database:
    
    def __init__(self, reset=False):
        if (reset == True):
            if os.path.exists("molecules.db"):
                os.remove("molecules.db")
        self.conn = sqlite3.connect('molecules.db')
        
                
    def create_tables( self ):
        cursor = self.conn.cursor()
        
        Elements = '''
        CREATE TABLE IF NOT EXISTS Elements (
            ELEMENT_NO INTEGER NOT NULL,
            ELEMENT_CODE VARCHAR(3) PRIMARY KEY NOT NULL,
            ELEMENT_NAME VARCHAR(32) NOT NULL,
            COLOUR1 CHAR(6) NOT NULL,
            COLOUR2 CHAR(6) NOT NULL,
            COLOUR3 CHAR(6) NOT NULL,
            RADIUS DECIMAL(3) NOT NULL);
        '''
        
        Atoms = '''
        CREATE TABLE IF NOT EXISTS Atoms (
            ATOM_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
            ELEMENT_CODE VARCHAR(3) NOT NULL,
            X DECIMAL(7, 4) NOT NULL,
            Y DECIMAL(7, 4) NOT NULL,
            Z DECIMAL(7, 4) NOT NULL,
            FOREIGN KEY (ELEMENT_CODE) REFERENCES Elements(ELEMENT_CODE));
        '''
        
        Bonds = '''
        CREATE TABLE IF NOT EXISTS Bonds (
            BOND_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
            A1 INTEGER NOT NULL,
            A2 INTEGER NOT NULL,
            EPAIRS INTEGER NOT NULL);
        '''
        
        Molecules = '''
        CREATE TABLE IF NOT EXISTS Molecules (
            MOLECULE_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
            NAME TEXT UNIQUE NOT NULL);
        '''
        
        MoleculeAtom = '''
        CREATE TABLE IF NOT EXISTS MoleculeAtom (
            MOLECULE_ID INTEGER NOT NULL,
            ATOM_ID INTEGER NOT NULL,
            PRIMARY KEY (MOLECULE_ID, ATOM_ID),
            FOREIGN KEY (MOLECULE_ID) REFERENCES Molecules(MOLECULE_ID),
            FOREIGN KEY (ATOM_ID) REFERENCES Atoms(ATOM_ID));
        '''
        
        MoleculeBond = '''
        CREATE TABLE IF NOT EXISTS MoleculeBond (
            MOLECULE_ID INTEGER NOT NULL,
            BOND_ID INTEGER NOT NULL,
            PRIMARY KEY (MOLECULE_ID, BOND_ID),
            FOREIGN KEY (MOLECULE_ID) REFERENCES Molecules(MOLECULE_ID),
            FOREIGN KEY (BOND_ID) REFERENCES Bonds(BOND_ID));
        '''
        
        cursor.execute(Elements)
        cursor.execute(Atoms)
        cursor.execute(Bonds)
        cursor.execute(Molecules)
        cursor.execute(MoleculeAtom)
        cursor.execute(MoleculeBond)
        
        self.conn.commit()
        
    def __setitem__( self, table, values ):
        cursor = self.conn.cursor()
        
        insertInto = f"INSERT INTO {table} VALUES {values}"
        
        cursor.execute(insertInto)
        
        self.conn.commit()
        
    def add_atom( self, molname, atom):
        cursor = self.conn.cursor()

        atomInfo = (atom.element, atom.x, atom.y, atom.z)
        
        #insertInto = f"INSERT INTO {Atoms} VALUES {atomInfo}"
        self.__setitem__("Atoms", atomInfo)
        #cursor.execute(insertInto)
        self.conn.commit()

        atomId = cursor.lastrowid
        cursor.execute(" SELECT MOLECULE_ID FROM Molecules WHERE name=?", (molname,))
        moleculeID = cursor.fetchone()[0]

        cursor.execute('''INSERT INTO MoleculeAtom (MOLECULE_ID, ATOM_ID) VALUES (?, ?),''', (moleculeID, atomId))
        self.conn.commit()

        
    def add_bond( self, molname, bond):
        cursor = self.conn.cursor()
        bondInfo = (bond.c_bond.a1, bond.c_bond.a2, bond.c_bond.epairs)
        
        self.__setitem__("Bonds", bondInfo)
        #cursor.execute(insertInto)
        self.conn.commit()
        
        bondId = cursor.lastrowid
        
        cursor.execute(" SELECT MOLECULE_ID FROM Molecules WHERE name=?", (molname,))
        moleculeID = cursor.fetchone()[0]
        
        cursor.execute('''INSERT INTO MoleculeBond (MOLECULE_ID, BOND_ID) VALUES (?, ?),''', (moleculeID, bondId))
        
        self.conn.commit()
        
    def add_molecule( self, name, fp ):    
        mol = MolDisplay.Molecule()
        mol.parse(fp)
        
        self.__setitem__("Molecules", name)
        
        numAtoms = 0
        numBonds = 0
    
        while (mol.get_atom(numAtoms)):
            numAtoms += 1
            
        while (mol.get_bond(numBonds)):
            numBonds += 1
            
        for i in range(numAtoms):
            self.add_atom(name, mol.get_atom(i))
        
        for i in range(numBonds):
            self.add_bond(name, mol.get_bond(i))
            
        
if __name__ == "__main__":
    db = Database(reset=True);
    db.create_tables();
    db['Elements'] = ( 'H', 1, 'Hydrogen', 'FFFFFF', '050505', '020202', 25 );
    db['Elements'] = ( 'C', 6, 'Carbon', '808080', '010101', '000000', 40 );
    db['Elements'] = ( 'N', 7, 'Nitrogen', '0000FF', '000005', '000002', 40 );
    db['Elements'] = ( 'O', 8, 'Oxygen', 'FF0000', '050000', '020000', 40 );
    fp = open( 'water-3D-structure-CT1000292221.sdf' );
    db.add_molecule( 'Water', fp );
    fp = open( 'caffeine-3D-structure-CT1001987571.sdf' );
    db.add_molecule( 'Caffeine', fp );
    fp = open( 'CID_31260.sdf' );
    db.add_molecule( 'Isopentanol', fp );
    print( db.conn.execute( "SELECT * FROM Elements;" ).fetchall() );
    print( db.conn.execute( "SELECT * FROM Molecules;" ).fetchall() );
    print( db.conn.execute( "SELECT * FROM Atoms;" ).fetchall() );
    print( db.conn.execute( "SELECT * FROM Bonds;" ).fetchall() );
    print( db.conn.execute( "SELECT * FROM MoleculeAtom;" ).fetchall() );
    print( db.conn.execute( "SELECT * FROM MoleculeBond;" ).fetchall() );
