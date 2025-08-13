# include "mol.h"

void atomset(atom * atom, char element[3], double * x, double * y, double * z)  {
    
    /* 
        Set the atom structures members with the corresponding parameters. the position in space of the atom is given by x, y, and z
        these will be copied into atom structures appropriate members. The name of the atom will also be copied into the appropraite 
        member.
    */

    atom->x = *x;
    atom->y = *y;
    atom->z = *z;
    //memcpy(atom->element, element, sizeof(char) * 3);
    strncpy(atom->element, element, 3);
}

void atomget(atom * atom, char element[3], double * x, double * y, double * z)  {
    
    /* 
        Set the fucntion parameters with atom structures the corresponding members. the position in space of the atom is given by 
        atom->x, atom->y, and atom->z these will be copied into the appropriate function parameters. The name of the atom given by
        atom->element will also be copied into the element parameter.
    */

    *x = atom->x;
    *y = atom->y;
    *z = atom->z;
    strncpy(element, atom->element, 3);
}

void bondset( bond *bond, unsigned short *a1, unsigned short *a2, atom **atoms, unsigned char *epairs ) {
    
    /*
        Set the bond structures members with the corresponding parameters. the atoms array will be refrenced by the two unsidned shorts
        a1, and a2 which will be copied into bond structures appropriate members. The type of bond will be spcified by epairs (ex single
        , double, or triple bonds) and stored in the corresponding bond structures member. The array's address will also be copied.
    */
    
    bond->a1 = *a1;
    bond->a2 = *a2;
    bond->epairs = *epairs;
    bond->atoms = *atoms;
    compute_coords(bond);
}

void bondget( bond *bond, unsigned short *a1, unsigned short *a2, atom **atoms, unsigned char *epairs ) {
    
    /* 
        Set the fucntion parameters with bond structures corresponding members. the two atom short are copied into place which is
        given by bond->a1, and bond->a2 these will be copied into the appropriate function parameters. The number of epairs will also
        be copied from the member into the fucntiuon parameter. Along with the correct address for atoms.
    */

    *a1 = bond->a1;
    *a2 = bond->a2;
    *epairs = bond->epairs;
    *atoms = bond->atoms;
}


void compute_coords( bond *bond )   {

    /*
        This function will compute the coordinats of the bond properly on the grid.
    */
    
    atom * a1 = &(bond->atoms[bond->a1]);
    atom * a2 = &(bond->atoms[bond->a2]);
    
    // Get the correct z-value which is will be the addition of the two bonds z and divide by 2.

    double divide = (double)2;
    double zVal = (a1->z + a2->z) / divide;
    
    double x1Val = a1->x;
    double x2Val = a2->x;
    double y1Val = a1->y;
    double y2Val = a2->y;

    double xDiff = x2Val - x1Val;
    double yDiff = y2Val - y1Val;

    double length = sqrt(pow(xDiff, 2) + pow(yDiff, 2));

    // Save the values in the members of the bond struct.

    bond->z = zVal;
    bond->x1 = x1Val;
    bond->x2 = x2Val;
    bond->y1 = y1Val;
    bond->y2 = y2Val;
    bond->len = length;
    bond->dx = xDiff / length;
    bond->dy = yDiff / length;
}

molecule *molmalloc( unsigned short atom_max, unsigned short bond_max ) {
    
    /*
        The goal of this function is to return a molecule structure's address, after it has been given the correct amount of space all
        given data will be intialized into the members and space will be malloced for all pointer arrays and chencked if it was properly
        done.
    */
    
    molecule * newMol = malloc(sizeof(molecule) * 1);

    // Check if malloc was a success.

    if (newMol == NULL) {
        return NULL;
    }

    // Access non pointer members and intialize them.

    newMol->atom_max = atom_max;
    newMol->atom_no = 0;

    // Allocate atom array members.

    newMol->atoms = malloc(sizeof(atom) * atom_max);
    newMol->atom_ptrs = (atom **)malloc(sizeof(atom *) * atom_max);

    // Check if malloc was a success for both pointer arrays.

    if (newMol->atoms == NULL)  {
        return NULL;
    }
    if (newMol->atom_ptrs == NULL)  {
        return NULL;
    }
    
    // Access non pointer members and intialize them.

    newMol->bond_max = bond_max;
    newMol->bond_no = 0;

    // Allocate bond array members.

    newMol->bonds = malloc(sizeof(bond) * bond_max);
    newMol->bond_ptrs = (bond **)malloc(sizeof(bond *) * bond_max);

    // Check if malloc was a success for both pointer arrays.

    if (newMol->bonds == NULL)  {
        return NULL;
    }
    if (newMol->bond_ptrs == NULL)  {
        return NULL;
    }

    return newMol;
}

molecule *molcopy( molecule *src )  {
    
    /*
        Copy the molecule given from src, all members must be copied along with all the data pointed to by src members into the 
        copiedMol.
    */
    
    molecule * copiedMol; 

    // If source is NULL return NULL, check to make sure you can even copy the source.

    if (src == NULL)    {
        return NULL;
    }

    // Allocate the molecule structure for the copied molecule along with all the members that can be changed.

    copiedMol = molmalloc(src->atom_max, src->bond_max);   
    
    // If molmalloc returns NULL we need to make sure the molcopy funciton also returns NULL before accessing the array pointer members.

    if (copiedMol == NULL)  {
        return NULL;
    }

    // Append the atoms and bonds to the molecule structure this happens after all malloc checks have passed.

    for (int i = 0; i < src->atom_no; ++i)  {
        molappend_atom(copiedMol, &(src->atoms[i]));
    }
        
    
    for (int i = 0; i < src->bond_no; ++i)  {
        molappend_bond(copiedMol, &(src->bonds[i]));
    }
        
    //copiedMol->atom_no = src->atom_no;
    //copiedMol->bond_no = src->bond_no;
  
    return copiedMol;
}

void molfree( molecule *ptr )   {
 
    /*
        Free all the pointer members and the pointer molecule structure itself. Also for good practive set the freed pointers to NULL.
    */
    
    // Free atoms array and set the pointer to NULL.

    free(ptr->atoms);
    ptr->atoms = NULL;

    // Free atom_ptrs array and set the pointer to NULL.

    free(ptr->atom_ptrs);
    ptr->atom_ptrs = NULL;
    
    // Free bonds array and set the pointer to NULL.

    free(ptr->bonds);
    ptr->bonds = NULL;
    
    // Free bond_ptrs array and set the pointer to NULL.

    free(ptr->bond_ptrs);
    ptr->bond_ptrs = NULL;

    // Free the struct itself and set it to NULL. 

    free(ptr);
    ptr = NULL;
}

void molappend_atom( molecule *molecule, atom *atom )   {
    
    /*
        The goal of this function is to populate the molecule structure with the information about the atoms in the molecule structure.
    */

    // Check if atom_no is equal to atom_max.
    
    if (molecule->atom_no == molecule->atom_max)    {
        
        // Increment atom_max accordingly.
        
        if (molecule->atom_max == 0)    {
            molecule->atom_max = 1;
        } 
        else    {
            molecule->atom_max *= 2;
        }
        
        // Allocate atoms and atom_ptrs arrays to makesure there is enough space for all the atoms including the new one.

        molecule->atoms = (struct atom *)realloc(molecule->atoms, sizeof(struct atom) * molecule->atom_max);
        molecule->atom_ptrs = (struct atom **)realloc(molecule->atom_ptrs, sizeof(struct atom *) * molecule->atom_max);

        /*
            Update atom_ptrs to point to the corresponding atoms in the new atoms array, to make sure there are no "misses" in the 
            atom_ptrs pointers.
        */ 

        for (int i = 0; i < molecule->atom_no; i++) {
            molecule->atom_ptrs[i] = &molecule->atoms[i];
        }
    }
    
    // Copy the atom to the first "empty" atom in atoms in the molecule. Followed by updating the atom_ptrs array.

    molecule->atoms[molecule->atom_no] = *atom;
    molecule->atom_ptrs[molecule->atom_no] = &molecule->atoms[molecule->atom_no];
    
    // Increment atom_no

    molecule->atom_no++;
}

void molappend_bond( molecule *molecule, bond *bond )   {
    
    /*
        The goal of this function is to populate the molecule structure with the information about the bonds in the molecule structure.
    */

    // Check if bond_no is equal to bond_max
    
    if (molecule->bond_no == molecule->bond_max)    {
        
        // Increment bond_max accordingly.

        if (molecule->bond_max == 0)    {
            molecule->bond_max = 1;
        }
        else    {
            molecule->bond_max *= 2;
        }
        
        // Allocate bonds and bond_ptrs arrays to makesure there is enough space for all the atoms including the new one.
        
        molecule->bonds = (struct bond *)realloc(molecule->bonds, sizeof(struct bond) * molecule->bond_max);
        molecule->bond_ptrs = (struct bond **)realloc(molecule->bond_ptrs, sizeof(struct bond *) * molecule->bond_max);

        /*
            Update bond_ptrs to point to the corresponding bonds in the new bonds array, to make sure there are no "misses" in the 
            bond_ptrs pointers.
        */ 

        for (int i = 0; i < molecule->bond_no; i++) {
            molecule->bond_ptrs[i] = &molecule->bonds[i];
        }
    }
    
    
    // Copy data from bond to the first "empty" bond in bonds in the molecule. Followed by updating the atom_ptrs array.

    molecule->bonds[molecule->bond_no] = *bond;
    molecule->bond_ptrs[molecule->bond_no] = &molecule->bonds[molecule->bond_no];
    
    // Increment atom_no

    molecule->bond_no++;
}

void molsort( molecule *molecule )  {

    /*
        The mol sort function will sor the atom_ptrs and bond_ptrs arrays in the molecule structure accordingly.
    */
    
    int numOfElements = (int)molecule->atom_no;
    int numOfBonds = (int)molecule->bond_no;
    
    // If statements cause program to only run qsort only if there are even atoms and bonds to sort.

    if (numOfElements > 0)  {

        // The base is molecule-atom_ptrs, nel = atom_no, size_t = atom * not atom, and the function to compare is passed.

        qsort(molecule->atom_ptrs, numOfElements, sizeof(atom *), atomCmp);
    }
    if (numOfBonds > 0) {

        // The base is molecule-bond_ptrs, nel = bond_no, size_t = bond * not bond, and the function to compare is passed.

        qsort(molecule->bond_ptrs, numOfBonds, sizeof(bond *), bondCmp);
    }
}

int atomCmp(const void *a, const void *b) {

    /*
        This function will be used to compare atoms in atom_ptrs array. It will compare induvisual atoms and get their z poisiton values
        to make it so that it will return an appropriate value depending on what the differnce between atom a and b are.  
    */
    
    // Type cast and assign the given void pointers that point to atom **.

    atom *aPtr = *(atom **)a;
    atom *bPtr = *(atom **)b;

    // This return statement allows for each of the 3 cases required by qsorts parameter to be performed.

    return (aPtr->z > bPtr->z) - (aPtr->z < bPtr->z);
}

int bondCmp(const void *a, const void *b) {
    
    /*
        This function will be used to compare bonds in bond_ptrs array. It will compare induvisual bonds and get their two atoms 
        average z poisiton value to make it so that it will return an appropriate value depending on what the differnce between
        them is.  
    */

    // Type cast and assign the given void pointers that point to bond **. Take the average z value of the atoms shared by that bond after.

    bond *aPtr = *(bond **)a;
    bond *bPtr = *(bond **)b;

    // This return statement allows for each of the 3 cases required by qsorts parameter to be performed.

    return (aPtr->z > bPtr->z) - (aPtr->z < bPtr->z);
}

void xrotation( xform_matrix xform_matrix, unsigned short deg ) {   
    
    /*
        This funciton will set the values in the xform_matrix according to a roatation of degeress around the x-axis

        | 1   0        0      |
        | 0   cosine   -sine  |
        | 0   sine     cosine |
    */
    
    double radians = deg * (M_PI / 180);
    double cosine = cos(radians);
    double sine = sin(radians);

    xform_matrix[0][0] = 1;
    xform_matrix[0][1] = 0;
    xform_matrix[0][2] = 0;
    xform_matrix[1][0] = 0;
    xform_matrix[1][1] = cosine;
    xform_matrix[1][2] = -1*sine;
    xform_matrix[2][0] = 0;
    xform_matrix[2][1] = sine;
    xform_matrix[2][2] = cosine;  
}

void yrotation( xform_matrix xform_matrix, unsigned short deg ) {   
    
    /*
        This funciton will set the values in the xform_matrix according to a roatation of degeress around the y-axis

        | cosine   0   sine   |
        | 0        1   0      |
        | -sine    0   cosine |
    */

    double radians = deg * (M_PI / 180);
    double cosine = cos(radians);
    double sine = sin(radians);

    xform_matrix[0][0] = cosine;
    xform_matrix[0][1] = 0;
    xform_matrix[0][2] = sine;
    xform_matrix[1][0] = 0;
    xform_matrix[1][1] = 1;
    xform_matrix[1][2] = 0;
    xform_matrix[2][0] = -sine;
    xform_matrix[2][1] = 0;
    xform_matrix[2][2] = cosine;  
}

void zrotation( xform_matrix xform_matrix, unsigned short deg ) {   
    
    /*
        This funciton will set the values in the xform_matrix according to a roatation of degeress around the z-axis

        | cosine   -sine    0 |
        | sine     cosine   0 |
        | 0        0        1 |
    */
    
    double radians = deg * (M_PI / 180);
    double cosine = cos(radians);
    double sine = sin(radians);

    xform_matrix[0][0] = cosine;
    xform_matrix[0][1] = -1*sine;
    xform_matrix[0][2] = 0;
    xform_matrix[1][0] = sine;
    xform_matrix[1][1] = cosine;
    xform_matrix[1][2] = 0;
    xform_matrix[2][0] = 0;
    xform_matrix[2][1] = 0;
    xform_matrix[2][2] = 1;  
}

void mol_xform( molecule *molecule, xform_matrix matrix )   {

    /*
        This funciton will apply the transformation matrix to all the atoms in the molecule.
    */
    
    int i = 0;
    int atom_no = (int)molecule->atom_no;
    int bond_no = (int)molecule->bond_no;
    double xPos;
    double yPos;
    double zPos;

    // Performs the vector multiplication to the x, y, z coordinates.

    for (i = 0; i < atom_no; i++)   {
        xPos =  molecule->atom_ptrs[i]->x;
        yPos =  molecule->atom_ptrs[i]->y;
        zPos =  molecule->atom_ptrs[i]->z; 

        molecule->atom_ptrs[i]->x = xPos * matrix[0][0] + yPos * matrix[0][1] + zPos * matrix[0][2];
        molecule->atom_ptrs[i]->y = xPos * matrix[1][0] + yPos * matrix[1][1] + zPos * matrix[1][2];
        molecule->atom_ptrs[i]->z = xPos * matrix[2][0] + yPos * matrix[2][1] + zPos * matrix[2][2];
    }

    for (i = 0; i < bond_no; i++)   {
        compute_coords(&(molecule->bonds[i]));
    }
}