#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

typedef struct Atoms {
  
  int atomNumber;
  char *atomName;
  char *residueName;
  int residueNumber;
  double x, y, z;
  float bfactor; 

} atom;

atom *load_atoms_data(FILE *pdb);

float *get_box_vectors(FILE *pdb);

float *find_box_center(float *vectors);

bool is_inside(double *coordinates, float *center, float cutoff);

int main(int argc, char *argv[])
{
  char *fileName; 
  char *outputName;
  char *resname = NULL;
  float cutoff;
  int counts = 0;

  for (int i = 1; i < argc; i++)
  {
    if (argv[i][0] == '-')
    {
      if (strcmp(argv[i], "-help") == 0)
      {
        printf("\nUsage: ./water_filter -file <input.pdb> -cutoff <cutoff_distance> -resname <WAT|TIP3|HOH|...> -output <output.pdb>\n\n");
        printf("Options:\n");
        printf("  -file      Path to the input PDB file.\n");
        printf("  -cutoff    Cutoff distance (in Ångströms) from the center of the box.\n");
        printf("  -resname   Residue name of water (e.g., SOL, WAT, HOH, TIP3).\n");
        printf("  -output    Path to the output PDB file (filtered result).\n");
        printf("  -help      Display this help message.\n\n");
        printf("Description:\n");
        printf("  This program removes water molecules (O and associated Hs) within the given cutoff\n");
        printf("  distance from the center of the simulation box. It assumes water atoms are ordered\n");
        printf("  as O followed by two H atoms in the PDB file.\n\n");
        printf("Example:\n");
        printf("  ./water_filter -file input.pdb -cutoff 5.0 -resname TIP3 -output filtered.pdb\n\n");
        return 0;
      }

    if (strcmp(argv[i], "-file") == 0)
      {
        if (i + 1 < argc)
        {
         fileName = argv[i + 1];
        }
      }

      if (strcmp(argv[i], "-resname") == 0)
      {
        if (i + 1 < argc)
        {
          resname = argv[i + 1];
        }
      }

      if (strcmp(argv[i], "-cutoff") == 0 )
      {
        cutoff = atof(argv[i + 1]);
      }
      if (strcmp(argv[i] , "-output") == 0)
      {  
        outputName = argv[i + 1];
      }
    } 
  } 
   
  FILE *pdb = fopen(fileName, "r");

  atom *atomArray = load_atoms_data(pdb);
  
  while(atomArray[counts].atomName != NULL)
  {
    counts++;
  }
  
  float *boxVectors = get_box_vectors(pdb);  

  float *center = find_box_center(boxVectors);

  double coordinates[3];

  for (int i = 0; i < counts;)
  {
 
    coordinates[0] = atomArray[i].x;
    coordinates[1] = atomArray[i].y;
    coordinates[2] = atomArray[i].z;
 

    if (atomArray[i].atomName[0] == 'O' && strncmp(atomArray[i].residueName, resname, strlen(resname)) == 0)
    {
      if (is_inside(coordinates, center, cutoff))
      {
        
        //assuming the SOL atoms are in order such OW HW HW 
        //then we just have to remove the line conatining OW
        //and the next 2 lines.

        atomArray[i].atomName = NULL;        
        atomArray[i + 1].atomName = NULL;
        atomArray[i + 2].atomName = NULL;
          
      } 
      i += 3;
    }
    else
    {
      i++;
    }
  }   
  
  FILE *outputFile = fopen(outputName, "w");

  for (int i = 0; i < counts; i++)
  {
    if (atomArray[i].atomName != NULL)
    {

      fprintf(outputFile, "ATOM  %5d %-4s%-3s %1s%4d    %8.3lf%8.3lf%8.3lf%6.2f%6.2f\n", 
            atomArray[i].atomNumber, atomArray[i].atomName, atomArray[i].residueName, 
            "", atomArray[i].residueNumber, atomArray[i].x, atomArray[i].y, atomArray[i].z,
            1.0, atomArray[i].bfactor);      
    }
  }
  
  for (int i = 0; i < counts; i++)
  {
    free(atomArray[i].atomName);
    free(atomArray[i].residueName);
  }

  free(atomArray);
  free(boxVectors);
  free(center);
  fclose(pdb);
  fclose(outputFile);
  return 0;
}

atom *load_atoms_data(FILE *pdb)
{
   
  char line[120];
  int atomCount = 0;

  while (fgets(line, sizeof(line), pdb) != NULL) 
  {
    if (strncmp(line, "ATOM", 4) == 0)
    {
      atomCount++;
    }
  }

  rewind(pdb);

  atom *atomArray = (atom *) malloc(atomCount * sizeof(atom));
  
  for (int i = 0; i < atomCount; i++)
  {
    
    atomArray[i].atomName = malloc(5 * sizeof(char));
    atomArray[i].residueName = malloc(6 * sizeof(char));

  }  


  int index = 0;

  char *type = malloc(7 * sizeof(char));
  char *atomName = malloc(5 * sizeof(char));
  int atomNumber;
  char *residueName = malloc(6 * sizeof(char));
  int residueNumber;
  double x, y, z;
  float one, bfactor;

  while (fgets(line, sizeof(line), pdb))
  {

    if (strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0)
    {
      atom *currentAtom = &atomArray[index];
     
      sscanf(line, "%s %d %s %s %d %lf %lf %lf %f %f", type, &atomNumber, atomName, 
             residueName, &residueNumber, &x, &y, &z, &one, &bfactor);

      atomArray[index].atomNumber = atomNumber;
      strcpy(atomArray[index].atomName, atomName);
      strcpy(atomArray[index].residueName, residueName);
      atomArray[index].residueNumber = residueNumber;
      atomArray[index].x = x;
      atomArray[index].y = y;
      atomArray[index].z = z;
      atomArray[index].bfactor = bfactor;
       
      index++; 
    }
  }

  return atomArray;
}

float *get_box_vectors(FILE *pdb)
{
  rewind(pdb);

  char line[120];
  char *cryst = malloc(5 * sizeof(char));
  char *strings = malloc(10 * sizeof(char));
  float *vectors = malloc(3 * sizeof(float));
  float *angles = malloc(3 * sizeof(float)); 

  for(int i = 0; i < 5; i ++)
  {
    fgets(line, sizeof(line), pdb);
    
    
    if (strncmp(line, "CRYST", 5) == 0 )
    {
      

      sscanf(line, "%s %f %f %f %f %f %f %s %s %s", cryst, &vectors[0], &vectors[1], &vectors[2], 
             &angles[0], &angles[1], &angles[2], &strings[0], &strings[1], &strings[2]);
    }
  }

  free(cryst);
  free(strings);
  
  return vectors;
}

float *find_box_center(float *vectors)
{
  float *center = malloc(3 * sizeof(float));

  for (int i = 0; i < 3; i++)
  {
    center[i] = vectors[i] / 2;
  }

  return center;
}

bool is_inside(double *coordinates, float *center, float cutoff)
{

  float distance = sqrt(powf(coordinates[0] - center[0], 2) + powf(coordinates[1] - center[1], 2) + powf(coordinates[2] - center [2], 2));
  
  return distance < cutoff;

}
