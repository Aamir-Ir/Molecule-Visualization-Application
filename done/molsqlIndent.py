#SQL.py

import os
import sqlite3

class Database:

	def __init__(self, reset=False):
		if (reset == True):
			os.remove("molecules.db")
    
    
    
		self.dataBase = sqlite3.connect('molecules.db')   
		
	def create_tables( self ):
		add = 0
  
  
	
	
  
	
	def __setitem__( self, table, values ):
		add = 0

	def add_atom( self, molname, atom ):
		add = 0
		
	def add_bond( self, molname, bond ):
		add = 0

	def add_molecule( self, name, fp ):
		add = 0

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


''' Before change
-- SQLite

CREATE TABLE Elements (
  ELEMENT_NO INTEGER NOT NULL,
  ELEMENT_CODE VARCHAR(3) PRIMARY KEY NOT NULL,
  COLOUR1 CHAR(6) NOT NULL,
  COLOUR2 CHAR(6) NOT NULL,
  COLOUR3 CHAR(6) NOT NULL,
  RADIUS DECIMAL(3) NOT NULL
);

CREATE TABLE Atoms (
  ATOM_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  ELEMENT_CODE VARCHAR(3) NOT NULL REFERENCES Elements(ELEMENT_CODE),
  X DECIMAL(7, 4) NOT NULL,
  Y DECIMAL(7, 4) NOT NULL,
  Z DECIMAL(7, 4) NOT NULL 
);

CREATE TABLE Bonds (
  BOND_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  A1 INTEGER NOT NULL REFERENCES Atoms(ATOM_ID),
  A2 INTEGER NOT NULL REFERENCES Atoms(ATOM_ID),
  EPAIRS INTEGER NOT NULL
);

CREATE TABLE Molecules (
  MOLECULE_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  NAME TEXT UNIQUE NOT NULL
);

CREATE TABLE MoleculeAtom (
  MOLECULE_ID INTEGER NOT NULL REFERENCES Molecules(MOLECULE_ID),
  ATOM_ID INTEGER NOT NULL REFERENCES Atoms(ATOM_ID),
  PRIMARY KEY (MOLECULE_ID, ATOM_ID)
);

CREATE TABLE MoleculeBond (
  MOLECULE_ID INTEGER NOT NULL REFERENCES Molecules(MOLECULE_ID),
  BOND_ID INTEGER NOT NULL REFERENCES Bonds(BOND_ID),
  PRIMARY KEY (MOLECULE_ID, BOND_ID)
);
'''
''' After Change
-- SQLite

CREATE TABLE Elements (
  ELEMENT_NO INTEGER NOT NULL,
  ELEMENT_CODE VARCHAR(3) PRIMARY KEY NOT NULL,
  ELEMENT_NAME VARCHAR(32) NOT NULL,
  COLOUR1 CHAR(6) NOT NULL,
  COLOUR2 CHAR(6) NOT NULL,
  COLOUR3 CHAR(6) NOT NULL,
  RADIUS DECIMAL(3) NOT NULL
);

CREATE TABLE Atoms (
  ATOM_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  ELEMENT_CODE VARCHAR(3) NOT NULL REFERENCES Elements(ELEMENT_CODE),
  X DECIMAL(7, 4) NOT NULL,
  Y DECIMAL(7, 4) NOT NULL,
  Z DECIMAL(7, 4) NOT NULL 
);

CREATE TABLE Bonds (
  BOND_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  A1 INTEGER NOT NULL,
  A2 INTEGER NOT NULL,
  EPAIRS INTEGER NOT NULL
);

CREATE TABLE Molecules (
  MOLECULE_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  NAME TEXT UNIQUE NOT NULL
);

CREATE TABLE MoleculeAtom (
  MOLECULE_ID INTEGER NOT NULL REFERENCES Molecules(MOLECULE_ID),
  ATOM_ID INTEGER NOT NULL REFERENCES Atoms(ATOM_ID),
  PRIMARY KEY (MOLECULE_ID, ATOM_ID)
);

CREATE TABLE MoleculeBond (
  MOLECULE_ID INTEGER NOT NULL REFERENCES Molecules(MOLECULE_ID),
  BOND_ID INTEGER NOT NULL REFERENCES Bonds(BOND_ID),
  PRIMARY KEY (MOLECULE_ID, BOND_ID)
);
'''
# molecules.sql


''' Final Correct one
-- SQLite

CREATE TABLE Elements (
  ELEMENT_NO INTEGER NOT NULL,
  ELEMENT_CODE VARCHAR(3) PRIMARY KEY NOT NULL,
  ELEMENT_NAME VARCHAR(32) NOT NULL,
  COLOUR1 CHAR(6) NOT NULL,
  COLOUR2 CHAR(6) NOT NULL,
  COLOUR3 CHAR(6) NOT NULL,
  RADIUS DECIMAL(3) NOT NULL
);

CREATE TABLE Atoms (
  ATOM_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  ELEMENT_CODE VARCHAR(3) NOT NULL,
  X DECIMAL(7, 4) NOT NULL,
  Y DECIMAL(7, 4) NOT NULL,
  Z DECIMAL(7, 4) NOT NULL,
  FOREIGN KEY (ELEMENT_CODE) REFERENCES Elements(ELEMENT_CODE)
);

CREATE TABLE Bonds (
  BOND_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  A1 INTEGER NOT NULL,
  A2 INTEGER NOT NULL,
  EPAIRS INTEGER NOT NULL
);

CREATE TABLE Molecules (
  MOLECULE_ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  NAME TEXT UNIQUE NOT NULL
);

CREATE TABLE MoleculeAtom (
  MOLECULE_ID INTEGER NOT NULL,
  ATOM_ID INTEGER NOT NULL,
  PRIMARY KEY (MOLECULE_ID, ATOM_ID),
  FOREIGN KEY (MOLECULE_ID) REFERENCES Molecules(MOLECULE_ID),
  FOREIGN KEY (ATOM_ID) REFERENCES Atoms(ATOM_ID)
);

CREATE TABLE MoleculeBond (
  MOLECULE_ID INTEGER NOT NULL,
  BOND_ID INTEGER NOT NULL,
  PRIMARY KEY (MOLECULE_ID, BOND_ID),
  FOREIGN KEY (MOLECULE_ID) REFERENCES Molecules(MOLECULE_ID),
  FOREIGN KEY (BOND_ID) REFERENCES Bonds(BOND_ID)
);
'''