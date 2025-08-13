from molecule import molecule

radius = { 'H': 25,
           'C': 40,
           'O': 40,
           'N': 40,
};

element_name = { 'H': 'grey',
                 'C': 'black',
                 'O': 'red',
                 'N': 'blue',
};

header = """<svg version="1.1" width="1000" height="1000" 
 xmlns="http://www.w3.org/2000/svg">""";

footer = """</svg>""";

offsetx = 500;
offsety = 500;

# Above is everything assignment 2 sheet specified to add.

# Atom class below. This class is used to essentially create an Atom and all the information required for it inorder to create an svg.

class Atom:

	# __init__ constructor method used so that the minimum fields are specified.

	def __init__(self, c_atom):
		self.c_atom = c_atom
		self.z = c_atom.z

	# __str__ method used to test and make sure all the values in each member of the atom struct were properly held by the self variable.

	def __str__(self):
		return "Element: {0} x: {1} y: {2} z: {3}".format(self.c_atom.element, self.c_atom.x, self.c_atom.y, self.c_atom.z)

	# svg method used to create the svg instructions based on each piece of information.

	def svg(self):

		# Computing the coordinates of the atoms along with the information for colour and element name.

		x = (self.c_atom.x * 100.0) + offsetx
		y = (self.c_atom.y * 100.0) + offsety
		r = radius[self.c_atom.element]
		colour = element_name[self.c_atom.element]
		return ' <circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s"/>\n' % (x, y, r, colour)

# Bond class below This class is used to essentially create an Bond and all the information required for it inorder to create an svg.

class Bond:

	# __init__ constructor method used so that the minimum fields are specified.

	def __init__(self, c_bond):
		self.c_bond = c_bond
		self.z = c_bond.z

	# __str__ method used to test and make sure all the values in each member of the bond struct were properly held by the self variable.

	def __str__(self):
		return "Epairs: %f x1: %f x2: %f y1: %f y2: %f z: %f len: %f dx: %f dy: %f" % (self.c_bond.epairs, self.c_bond.x1, self.c_bond.x2, self.c_bond.y1, self.c_bond.y2, self.c_bond.z, self.c_bond.len, self.c_bond.dx, self.c_bond.dy)

	# svg method used to create the svg instructions based on each piece of information.

	def svg(self):

		# Assign values for readable code.

		dx = self.c_bond.dx
		dy = self.c_bond.dy
		x1 = self.c_bond.x1
		y1 = self.c_bond.y1
		x2 = self.c_bond.x2
		y2 = self.c_bond.y2

		'''
			Compute the Bonds. Essentially apply the grid by multiplying each value by 100 then the add/subtract (depends on which way you mover)
		 	the change (dx/dy) * grid then finally the offset. 
		'''

		x1_l1 = (x1 * 100) - (dy * 10) + offsety
		y1_l1 = (y1 * 100) + (dx * 10) + offsetx
		x1_h2 = (x1 * 100) + (dy * 10) + offsety
		y1_h2 = (y1 * 100) - (dx * 10) + offsety

		x2_l1 = (x2 * 100) - (dy * 10) + offsety
		y2_l1 = (y2 * 100) + (dx * 10) + offsetx
		x2_h2 = (x2 * 100) + (dy * 10) + offsety
		y2_h2 = (y2 * 100) - (dx * 10) + offsetx

		# Note the order you display the points do not matter so you are good.

		return ' <polygon points="%.2f,%.2f %.2f,%.2f %.2f,%.2f %.2f,%.2f" fill="green"/>\n' % (x1_l1, y1_l1, x1_h2, y1_h2, x2_h2, y2_h2, x2_l1, y2_l1)

'''
	Molecule subclass of molecule class below. This class is used to essentially create a Molecule and it does this by parsing a file for information about atoms
	and bonds and uses the classes above to create the SVG for them then finally creates the final instructions for svg to create the entire molecule on the 
	webserver.
'''

class Molecule(molecule):

	# __str__ method used to test all the output.

	def __str__(self):

		for i in range(numAtoms):
			atoms = self.get_atom(i)
			print("Atom %d: " % (i + 1) + atoms.element, atoms.x, atoms.y, atoms.z)

		for i in range(numBonds):
			bonds = self.get_bond(i)
			print("Bond %d: " % (i + 1) + str(bonds.a1), str(bonds.a2), str(bonds.epairs), str(bonds.x1), str(bonds.y1), str(bonds.x2), str(bonds.y2), str(bonds.len), str(bonds.dx), str(bonds.dy))

	'''
		The svg method is used to generate a string that will be encoded so that the webform will be able to display the molecule, it does this by calling the
		atom, and bond classes then finally collects the svg method returns from them in the increasing z value order. Finally there will be an output for
		the svg returned by a string containing the instructions to build the moelecule.
	'''

	def svg(self):

		# Start with the header. then collect the info on the a values in increasing order for the atoms and bonds in 2 lists.

		string = header

		self.sort()

		atomZ = []
		bondZ = []

		for x in range(numAtoms):
			atomZ.append(self.get_atom(x).z)

		for y in range(numBonds):
			bondZ.append(self.get_bond(y).z)

		#print(atomZ)
		#print(bondZ)

		# Traverse the two lists and compare the z values and in increasing order collect he svg's for the atom or bond in the correct order.

		i = 0
		j = 0

		while (i < len(atomZ) or j < len(bondZ)):

			# Break conditions for the lists as to avoid out of index error.

			if ( i == len(atomZ)):
				break
			if (j == len(bondZ)):
				break

			# Applying the algorithm spcified in assignment 2.

			a = Atom(self.get_atom(i))
			b = Bond(self.get_bond(j))

			if a.z < b.z:
				atomStr = a.svg()
				#print("i = " + str(i))
				#print(atomStr)
				string += atomStr
				i += 1
			else:
				bondStr = b.svg()
				#print("j = " + str(j))
				#print(bondStr)
				string += bondStr
				j += 1

		# Finally add the svg strings for the last remaining atoms or bonds.

		pos1 = i
		pos2 = j

		inting1 = len(atomZ) - pos1
		inting2 = len(bondZ) - pos2

		if (i == len(atomZ)):
			for pos2 in range(inting2):
				b = Bond(self.get_bond(j))
				bondStr = b.svg()
				j += j
				#print("j = " + str(j))
				string += bondStr
		if (j == len(bondZ)):
			for pos1 in range(inting1):
				a = Atom(self.get_atom(i))
				atomStr = a.svg()
				i += 1
				#print("i + " + str(i))
				string += atomStr

		string += footer

		#print(string)
		return string

	# The parse method is used to collect information from the input file.

	def parse(self, file):

		# Skip the first 3 lines.

		next(file)
		next(file)
		next(file)

		# Read the information of the molecule referring to how many atoms and bonds there are.

		atomBondInfo = file.readline().strip()

		# Set the values for the num of atoms and bonds and provide a scope for the atoms and bonds to help other classes use the num of atoms and bonds.

		global numAtoms
		numAtoms = int(atomBondInfo.split()[0])
		global numBonds
		numBonds = int(atomBondInfo.split()[1])

		# Collect all the data for atoms line by line, use unpacking to make it easier.

		for i in range(numAtoms):
			atomData = file.readline().strip()
			x, y, z, element = float(atomData.split()[0]), float(atomData.split()[1]), float(atomData.split()[2]), atomData.split()[3]
			#print(x, y, z, element)
			self.append_atom(element, x, y, z)

		# Collect all the data for bonds line by line, use unpacking to make it easier.

		for i in range(numBonds):
			bondData = file.readline().strip()
			a1, a2, bondType = int(bondData.split()[0]), int(bondData.split()[1]), int(bondData.split()[2])
			#print(str(a1), str(a2), str(bondType))
			self.append_bond(a1, a2, bondType)


