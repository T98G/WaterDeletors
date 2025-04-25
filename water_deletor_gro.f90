program water_deletor

  implicit none
  
  integer :: i, stats
  integer :: numAtoms, newNumAtoms
  integer, pointer :: residues(:)  
  real :: x, y, z, distance, radius
  real :: boxVectors(3), center(3)
  character(len=100) :: gro, outputFile
  logical :: help
  character(len=30), dimension(30) :: args
  character(len=80) :: line, header
  character(len=6) :: waterName

  type AtomData
    character(len=4) :: name
    character(len=6) :: residueName
    integer :: residueNumber
    integer :: number
    real :: position(3)
  end type AtomData

  type(AtomData), allocatable :: atomArray(:)
  type(AtomData), allocatable :: filtered(:)

  help = .false.

  do i = 1, 30
    
    call getarg(i, args(i))

  end do

  do i = 1, 30
    
    if (args(i) == "-gro") then
      
      gro = trim(args(i + 1))
    
    end if

    if (args(i) == "-wres") then

      waterName = trim(args(i + 1))


    end if

    if (args(i) == "-radius") then

      read(args(i + 1), *) radius

    end if

    if  (args(i) == "-o") then

      outputFile = trim(args(i + 1))

    end if

    if (args(i) == "-help") then

      help = .true.

    end if

  end do

  if (help) then

    write(*, *) "water_deletor - Removes water molecules within a given radius from the center of a simulation box."
    write(*, *)
    write(*, *) "Usage:"
    write(*, *) "  water_deletor -gro input.gro -radius R -o output.gro"
    write(*, *)
    write(*, *) "Options:"
    write(*, *) "  -gro <filename>      Input GROMACS .gro structure file."
    write(*, *) "  -radius <value>      Radius (in nm) around the center of the box to delete water molecules."
    write(*, *) "  -wres <string>       The residue name of the water molecules in the file"
    write(*, *) "  -o <filename>        Output .gro file with selected water molecules removed."
    write(*, *) "  --help               Display this help message and exit."
    write(*, *)
    write(*, *) "Example:"
    write(*, *) "  water_deletor -gro input.gro -radius 1.0 -o cleaned.gro"

    stop

  end if


  open(unit = 1, file = gro, status = "old", access = "sequential", form = "formatted", action = "read")

  !Read the header and atom number

  read(1, "(A)", iostat = stats) line
  
  header = line
  
  read(1, "(A)", iostat = stats) line
  
  read(line, *) numAtoms
    

  allocate(atomArray(numAtoms))

  !Not rewinding since I want to start reading from line 2

  do i = 1, numAtoms 

    read(1, "(A)", iostat = stats) line
    
    read(line(1:5), *) atomArray(i)%residueNumber
    atomArray(i)%residueName = trim(line(6:11))
    atomArray(i)%name = trim(line(12:15))
    read(line(16:20), *) atomArray(i)%number
    read(line(24:28), *) atomArray(i)%position(1)
    read(line(32:36), *) atomArray(i)%position(2)
    read(line(40:44), *) atomArray(i)%position(3)
 
  end do

  read(1, "(A)", iostat = stats) line
  read(line(4:10), *) boxVectors(1)
  read(line(14:20), *) boxVectors(2)
  read(line(24:30), *) boxVectors(3)
  
  !close gro file
  close(1)

  center = boxVectors/2

  allocate(residues(numAtoms))

  do i = 1, numAtoms
    
    x = atomArray(i)%position(1)
    y = atomArray(i)%position(2)
    z = atomArray(i)%position(3)

    distance = sqrt((x - center(1)) ** 2 + (y - center(2)) ** 2 + (z - center(3)) ** 2)

    if (distance < radius .and. atomArray(i)%residueName(1:4) == waterName) then

      residues(i) = atomArray(i)%residueNumber

    end if
  
  end do
  
  allocate(filtered(numAtoms))

  do i = 1, numAtoms
  

    if (.not. any(residues == atomArray(i)%residueNumber)) then
         
        filtered(i)%name = atomArray(i)%name
        filtered(i)%number = atomArray(i)%number
        filtered(i)%residueName = atomArray(i)%residueName
        filtered(i)%residueNumber = atomArray(i)%residueNumber
        filtered(i)%position(1) = atomArray(i)%position(1)
        filtered(i)%position(2) = atomArray(i)%position(2)
        filtered(i)%position(3) = atomArray(i)%position(3)
    end if
  end do

  open(unit = 2, file = outputFile, status = "replace", action = "write")
  
  newNumAtoms = count(filtered%number /= 0)


  write(2, "(A)") header
  write(2, "(I8)") newNumAtoms

  do i = 1, numAtoms

    if (filtered(i)%number /= 0) then

      write(2, "(I5, A5, A1, A4, I5, f8.3, f8.3, f8.3)") &
        filtered(i)%residueNumber, filtered(i)%residueName, " ", trim(filtered(i)%name), &
        filtered(i)%number, filtered(i)%position(1), filtered(i)%position(2), &
        filtered(i)%position(3)
    
    end if

  end do
  
  write(2, "(A6, 3(f8.3))") "     ", boxVectors(1), boxVectors(2), boxVectors(3)

  close(2)
  
  deallocate(atomArray)
  deallocate(residues)
  deallocate(filtered)

end program water_deletor 
