! 
!     Copyright (C) 2021  Chee Kwan Gan (ihpcganck@gmail.com)
!     Copyright (C) 2020  Chee Kwan Gan (ihpcganck@gmail.com)
! 
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
! 

module extra
  use commod
  implicit none
  integer,parameter :: inputu   =11
  integer,parameter :: NNu      =12
  integer,parameter :: vaspu    =13
  integer,parameter :: siestau  =14
  integer,parameter :: pwposcaru=15
  integer,parameter :: xyzu     =16
  integer,parameter :: fleuru   =17
end module extra

program sample
  use extra
  implicit none
  type(supercell) :: s,t,f
  character(len=DSL) :: filestem
  character(len=DSL) :: inputfile,outputfile,tmp
  character(len=DSL) :: inputformat,outputformat
  character(len=DSL) :: InputLengthUnit,OutputLengthUnit,ocf
  integer :: n1,n2,n3
  integer :: atmtype,zord(zordmax)
  integer :: nmax
  logical :: do_NN
  character(len=DSL) :: qtablepath

  call init_static_sg_database(qtablepath)

  call fargv(1,inputfile)
  if(trim(inputfile) == 'help') then
    write(*,'(A)') ' babel (1) input file (vasp,arc,car,STRUCT_OUT,gjf,fdf,pdb)'
    write(*,'(A)') '   (2) input unit (ang,bohr,...) (3) output file format (xsf,cif,vasp,fdf,gjf,gulp,pdb)'
    write(*,'(A)') '   (4) output unit (5-7) n1,n2,n3 (8) output coordinate format (frac/abs)'
    write(*,*) 'for example:'
    write(*,*) './babel cell.vasp bohr cif ang 1 2 1 frac'
    stop
  endif
  call fargn(8)
  call fargv(2,InputLengthUnit)
  call fargv(3,outputformat)
  call fargv(4,OutputLengthUnit)
  call fargv(5,n1)
  call fargv(6,n2)
  call fargv(7,n3)
  call fargv(8,ocf)

  tmp = 'echo: babel '//trim(inputfile)//' '//trim(InputLengthUnit)
  tmp = trim(tmp)//' '//trim(outputformat)//' '//trim(OutputLengthUnit)
  tmp = trim(tmp)//' '//trim(N2str(n1))//' '//trim(N2Str(n2))//' '//trim(N2str(n3))//' '//trim(ocf)
  write(*,'(A)') trim(tmp)

  call read_struc(inputu,inputfile,filestem,inputformat,s)

  if(trim(InputLengthUnit) == trim(OutputLengthUnit) .and. trim(inputformat) == trim(outputformat)) then
    write(*,*) 'Warning: same InputLengthUnit and OutputLengthUnit. same inputformat and outputformat'
    outputfile = 'reformat'//'.'//trim(outputformat)
  else
    outputfile = trim(filestem)//'.'//trim(outputformat)
  endif

  if(trim(outputformat) == 'pwposcar') then
    write(*,'(A)') 'outputfiles are '//trim(outputfile)//'1'//' and '//trim(outputfile)//'2'
  else
    write(*,'(A)') 'babel: outputfile = '//trim(outputfile)
  endif

  if(trim(inputfile) == trim(outputfile)) then
    write(*,*) 'warning: same file names for inputfile and outputfile'
    outputfile = 'new-'//trim(outputfile)
    write(*,*) 'Since you are dealing the same format, you have to use a new filename.'
    write(*,*) 'The new file name is ',trim(outputfile)
  endif

  do_NN = .false.
  if(do_NN) then
    write(*,*) 'Table of nearest neighbours ...'
    open(unit=NNu,file='NN-dist-tab.dat',status='replace')
    nmax=5000
    call nearest_neighbor_tab(NNu,s,2,2,2,nmax)
    close(NNu)
  endif

  if(trim(InputLengthUnit)=='bohr' .and. trim(OutputLengthUnit) == 'ang') then
    call cp_cell_with_scaling(s,t,Bohr2Ang)
  else if(trim(InputLengthUnit)=='ang' .and. trim(OutputLengthUnit) == 'bohr') then
    call cp_cell_with_scaling(s,t,Ang2Bohr)
  else if(trim(InputLengthUnit)=='ang' .and. trim(OutputLengthUnit) == 'ang') then
    call cp_cell_with_scaling(s,t,1d0)
  else if(trim(InputLengthUnit)=='bohr' .and. trim(OutputLengthUnit) == 'bohr') then
    call cp_cell_with_scaling(s,t,1d0)
  else
    stop 'err: check your units under inputformat == vasp'
  endif

  call get_nsp_zord(t,atmtype,zord(1))

  call enlarge_supercell(t,n1,n2,n3,f)
  write(*,*) 'density is ',crys_density(t),' g/cm^3'

  f%fracabs = trim(ocf)
  f%sg_label = trim(s%sg_label)
  if(trim(outputformat)=='xsf') then
    call supercell_2_xcrysden(trim(outputfile),f)
    write(*,*) 'xcrysden --xsf '//trim(outputfile)
  else if(trim(outputformat) == 'cif') then
    call supercell_2_cif(trim(outputfile),atmtype,zord(1),f)
  else if(trim(outputformat) == 'vasp') then
    call supercell_2_vasp(vaspu,trim(outputfile),atmtype,zord(1),f)
  else if(trim(outputformat) == 'fdf') then
    call supercell_2_siesta_fdf(siestau,trim(outputfile),atmtype,zord(1),f)
  else if(trim(outputformat) == 'pwposcar') then
    call supercell_2_pwposcar(pwposcaru,trim(outputfile),atmtype,zord(1),f)

  else if(trim(outputformat) == 'gulp') then
    call supercell_2_gulp(trim(outputfile),atmtype,zord(1),f)
  else if(trim(outputformat) == 'xyz') then
    call supercell_2_xyz(xyzu,trim(outputfile),f)
  else if(trim(outputformat) == 'fleur') then
    call supercell_2_fleur(fleuru,trim(outputfile),atmtype,zord(1),f)
  else
    write(*,*) 'supported formats: xsf,cif,vasp,pwposcar,fdf,gjf,gulp'
    stop 'err: only a few format are supported.'
  endif
  write(*,'(A)') 'babel: '//trim(outputfile)//' is generated'

end program sample
