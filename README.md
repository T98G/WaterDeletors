# WaterDeletors

# waterdeletor.c

This program is designed to remove water molecules (oxygen and associated hydrogens) from a PDB file that are within a specified cutoff distance from the center of the simulation box. It is particularly useful for processing water molecules that are near the box center, which may be unnecessary for certain types of simulations or analyses.

## Features

- Removes water molecules (O and its associated H atoms) within a given cutoff distance from the center of the simulation box.
- Supports specifying the residue name (e.g., TIP3, WAT, HOH) for water molecules.
- Filters atoms from a PDB file based on the defined criteria.
- Outputs the filtered PDB file to a user-specified location.

## Usage

### Command Line Syntax

```
./water_filter -file <input.pdb> -cutoff <cutoff_distance> -resname <WAT|TIP3|HOH|...> -output <output.pdb>
```

### Options

- `-file <input.pdb>`: Path to the input PDB file containing the structure.
- `-cutoff <cutoff_distance>`: Cutoff distance (in Ångströms) from the center of the box to filter water molecules.
- `-resname <residue_name>`: Residue name of the water molecules (e.g., SOL, WAT, HOH, TIP3).
- `-output <output.pdb>`: Path to the output PDB file containing the filtered result.
- `-help`: Displays a help message with usage instructions.

### Example

To run the program with an input PDB file `input.pdb`, a cutoff distance of 5.0 Å, filtering water molecules with residue name `TIP3`, and output the filtered structure to `filtered.pdb`:

```
./water_filter -file input.pdb -cutoff 5.0 -resname TIP3 -output filtered.pdb
```

## Description

The program reads an input PDB file, extracts atom information, and checks if water molecules (with the specified residue name) are within the cutoff distance from the center of the simulation box. If a water molecule is within this distance, it removes the oxygen atom and its two associated hydrogen atoms from the structure.

### Key Functions:

- `load_atoms_data(FILE *pdb)`: Loads the atom data from the PDB file.
- `get_box_vectors(FILE *pdb)`: Retrieves the vectors representing the dimensions of the simulation box.
- `find_box_center(float *vectors)`: Calculates the center of the simulation box based on the box vectors.
- `is_inside(double *coordinates, float *center, float cutoff)`: Determines if an atom is within the cutoff distance from the center of the box.

## Memory Management

The program allocates memory dynamically for the atom array and various other data structures. It is important to free this memory before the program exits.

## Compilation

You can compile the program using a C compiler, such as `gcc`:

```bash
gcc -o water_filter water_filter.c -lm
```

Ensure that the math library (`-lm`) is linked during compilation for the proper functioning of mathematical operations like `sqrt()`.

## License

This program is distributed under the MIT License. See the LICENSE file for more details.

---

# WaterDeletor_gro.f90

The **Water Deletor** program removes water molecules (specified by residue name) within a defined radius from the center of the simulation box in a GROMACS `.gro` file. This tool is useful for cleaning up water molecules that are too close to the box center, which may be unnecessary for certain simulations or analyses.

## Features

- Reads GROMACS `.gro` files.
- Removes water molecules within a given radius from the center of the simulation box.
- Supports specifying the water residue name (e.g., `WAT`, `TIP3`).
- Outputs a new `.gro` file with the filtered result.

## Usage

### Command Line Syntax

```
water_deletor -gro <input.gro> -radius <radius> -wres <water_residue> -o <output.gro>
```

### Options

- `-gro <filename>`: Path to the input GROMACS `.gro` structure file.
- `-radius <value>`: Radius (in nm) around the center of the box within which water molecules will be deleted.
- `-wres <water_residue>`: The residue name of the water molecules (e.g., `WAT`, `TIP3`).
- `-o <filename>`: Path to the output `.gro` file with selected water molecules removed.
- `-help`: Displays a help message with usage instructions.

### Example

To run the program with an input GROMACS `.gro` file `input.gro`, a radius of 1.0 nm, removing water molecules with residue name `WAT`, and output the cleaned structure to `cleaned.gro`:

```
water_deletor -gro input.gro -radius 1.0 -wres WAT -o cleaned.gro
```

## Description

The program reads an input GROMACS `.gro` file, extracts the atom data, and calculates the distance of each water molecule from the center of the simulation box. If a water molecule is within the specified radius, it is removed from the output structure. The remaining atoms are written to a new `.gro` file.

### Key Functions:

- **Reading the `.gro` file**: The program reads atom data, including the residue name, atom name, atom number, and coordinates.
- **Filtering water molecules**: The program checks if each water molecule is within the specified radius from the box center and removes it if it is.
- **Writing the filtered output**: The program writes the remaining atoms (excluding the filtered water molecules) to a new `.gro` file.

### Program Logic:

1. The program first parses the command-line arguments to extract the input file, water residue name, radius, and output file name.
2. It then reads the `.gro` file, extracting the atom data and box dimensions.
3. For each atom, it calculates the distance to the center of the box and marks the water molecules within the specified radius for removal.
4. The filtered atom list is written to the output `.gro` file.

## Compilation

You can compile the program using a Fortran compiler, such as `gfortran`:

```bash
gfortran -o water_deletor water_deletor.f90
```

## License

This program is distributed under the MIT License. See the LICENSE file for more details.

---

