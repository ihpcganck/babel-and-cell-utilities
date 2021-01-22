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

module commod
  use sgsym
  use derived_constants
  implicit none

  interface fargv
    module procedure getarg_str,getarg_int,getarg_double
  end interface

contains



  subroutine getarg_str(n,str)
    integer :: n
    character(len=*) :: str
    call get_command_argument(n,str)
  end subroutine getarg_str

  subroutine getarg_int(n,k)
    integer :: n, k
    character(len=100) :: tmpstr
    call get_command_argument(n,tmpstr)
    read(unit=tmpstr,FMT='(I10)') k
  end subroutine getarg_int

  subroutine getarg_double(n,k)
    integer :: n
    real(double) :: k
    character(len=100) :: tmpstr
    call get_command_argument(n,tmpstr)
    read(unit=tmpstr,FMT='(D22.15)') k
  end subroutine getarg_double

  function str2N(str) result(n)
    character(len=*) :: str
    integer :: n
    read(str,'(I10)') n
  end function str2N

  function N2str(i) result (str)
    integer :: i,j,tmpv,intarr(20),nd,digit
    character(len=20) :: str

    if(i < 0) then
      tmpv = -i
    else
      tmpv = I
    endif

    nd = 0
    do
      digit = modulo(tmpv,10)
      nd = nd + 1
      intarr(nd) = digit
      tmpv = tmpv/10
      if(tmpv == 0) exit
    enddo
    str = ' '
    do j = 1, nd
      str(j:j) = digitarr(intarr(nd-j+1))
    enddo
    if(i < 0) then
      nd = nd + 1
      do j = nd,2,-1
        str(j:j) = str(j-1:j-1)
      enddo
      str(1:1) = '-'
    endif
  end function N2str

  subroutine fargn(n)
    integer :: n
    integer :: m

    m = command_argument_count()
    if(n /= m) then
      write(*,'(A,I4,A)') 'we need ',n,' arguments'
      write(*,'(A)') 'err: check number of arguments'
      stop 1
    endif
  end subroutine fargn

  function num_of_strings(tmp)
    character(len=*) :: tmp
    integer :: num_of_strings
    character(len=DSL) :: line
    integer,parameter :: max_num_cols = 100

    character(len=100) :: test_array(max_num_cols)
    integer :: i,io

    line = trim(tmp)
    write(*,'(A)') 'line is '//trim(line)
    do i = 1, max_num_cols

      read(line,*,iostat=io) test_array(1:i)

      if(io /= 0) then
        exit
      endif
    enddo
    num_of_strings = i-1
    if(num_of_strings == max_num_cols) then

    endif
  end function num_of_strings

  function num_of_integers(tmp)
    character(len=*) :: tmp
    integer :: num_of_integers
    character(len=DSL) :: line
    integer,parameter :: max_num_cols = 100
    integer :: test_array(max_num_cols)
    integer :: i,io

    line = trim(tmp)
    write(*,'(A)') 'line is '//trim(line)
    do i = 1, max_num_cols

      read(line,*,iostat=io) test_array(1:i)
      if(io /= 0) then
        exit
      endif
    enddo
    num_of_integers = i-1
    if(num_of_integers == max_num_cols) then
      write(*,*) 'count_int_number= ',num_of_integers
      write(*,*) 'err: we may need to increase max_num_cols.'
      stop 1
    endif
  end function num_of_integers

  subroutine getlenang(a,L)
    real(double) :: a(3,3),L(6)
    real(double) :: sintheta,costheta

    L(1) = vecmag3(a(1:3,1))
    L(2) = vecmag3(a(1:3,2))
    L(3) = vecmag3(a(1:3,3))

    costheta = dotprod3(a(1:3,2),a(1:3,3))/(L(2)*L(3))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(4) = atan2(sintheta,costheta)*rad2deg

    costheta = dotprod3(a(1:3,3),a(1:3,1))/(L(3)*L(1))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(5) = atan2(sintheta,costheta)*rad2deg

    costheta = dotprod3(a(1:3,1),a(1:3,2))/(L(1)*L(2))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(6) = atan2(sintheta,costheta)*rad2deg
  end subroutine getlenang

  subroutine real2recip(reallatt,reciplatt)
    real(double) :: reallatt(3,3),reciplatt(3,3)
    real(double) :: vol,tmp2(3,3)

    vol = det3(reallatt(1,1))
    if(vol < 0) then
      write(*,*) 'vol = ',vol
      write(*,'(A)') 'err: vol is negative. Check reallatt.'
      stop 1
    endif

    tmp2 = inv3x3(reallatt(1,1))

    reciplatt(1:3,1) = tmp2(1,1:3)
    reciplatt(1:3,2) = tmp2(2,1:3)
    reciplatt(1:3,3) = tmp2(3,1:3)
  end subroutine real2recip

  function dotprod3(A,B)
    real(double) :: A(1:3),B(1:3),dotprod3
    dotprod3 = dotprodn(3,a(1),b(1))
  end function dotprod3

  function dotprodn(n,A,B)
    integer :: n,i
    real(double) :: A(n),B(n),dotprodn
    dotprodn = zero
    do i = 1, n
      dotprodn = dotprodn + A(i)*B(i)
    enddo
  end function dotprodn

  function tripprod(A)
    real(double) :: A(3,3),tripprod
    tripprod = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - &
               A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + &
               A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  end function tripprod

  function det3(a)
    real(double) :: det3,a(3,3)
    det3 = tripprod(a)
  end function det3

  subroutine frac2abs(reallatt,FracC,AbsC)
    real(double) :: FracC(1:3),AbsC(1:3),reallatt(3,3)
    integer :: i,j
    do i = 1, 3
      AbsC(i) = 0
      do j = 1, 3
        AbsC(i) = AbsC(i) + reallatt(i,j)*FracC(j)
      enddo
    enddo
  end subroutine frac2abs

  subroutine abs2frac(reciplatt,AbsC,FracC)
    real(double) :: AbsC(3),FracC(3),reciplatt(3,3)
    integer :: i
    do i = 1, 3
      FracC(i) = dotprod3(AbsC(1:3),reciplatt(1:3,i))
    enddo
  end subroutine abs2frac

  subroutine enlarge_supercell(s1,n1,n2,n3,s2)
    type(supercell) :: s1,s2
    integer :: n1,n2,n3,i,j,k,m,ind
    real(double) :: newp(3),x(3)

    s2%sg_label = trim(s1%sg_label)
    s2%a(1:3,1) = n1*s1%a(1:3,1)
    s2%a(1:3,2) = n2*s1%a(1:3,2)
    s2%a(1:3,3) = n3*s1%a(1:3,3)
    call getlenang(s2%a,s2%la)
    call real2recip(s2%a,s2%b)
    s2%n = s1%n*n1*n2*n3
    ind = 0
    do i = 1, n1
      do j = 1, n2
        do k = 1, n3
          do m = 1, s1%n
            ind = ind + 1

            if(ind > maxnum) then
              write(*,*) 'In enlarge_supercell: '
              write(*,*) 'ind,maxnum = ',ind,maxnum
              write(*,'(A)') 'err: increase maxnum.'
              stop 1
            endif
            newp = (/i-1,j-1,k-1/) + s1%at(m)%f

            call frac2abs(s1%a,newp,x)
            call abs2frac(s2%b,x,s2%at(ind)%f)

            call frac2abs(s2%a,s2%at(ind)%f,s2%at(ind)%ac)
            s2%at(ind)%z = s1%at(m)%z
            s2%at(ind)%chg = s1%at(m)%chg
          enddo
        enddo
      enddo
    enddo
  end subroutine enlarge_supercell

  function vecmag3(x)
    real(double) :: vecmag3,x(3)
    vecmag3 = sqrt(dotprod3(x,x))
  end function vecmag3

  function inv3x3(a) result(inv)
    real(double) :: a(3,3),inv(3,3)
    real(double) :: det,absdet

    inv(1,1) = a(2,2)*a(3,3)-a(2,3)*a(3,2)
    inv(2,1) = -(a(2,1)*a(3,3) - a(2,3)*a(3,1))
    inv(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
    inv(1,2) = -(a(1,2)*a(3,3) - a(1,3)*a(3,2))
    inv(2,2) = a(1,1)*a(3,3) - a(1,3)*a(3,1)
    inv(3,2) = -(a(1,1)*a(3,2) - a(1,2)*a(3,1))
    inv(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
    inv(2,3) = -(a(1,1)*a(2,3) - a(1,3)*a(2,1))
    inv(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

    det = a(1,1)*inv(1,1)+a(1,2)*inv(2,1)+a(1,3)*inv(3,1)
    absdet = abs(det)

    if(absdet < 1.0d-13) then
      write(*,*) 'absdet = ',absdet
      write(*,'(A)') 'err: det is too small in inv3x3.'
      stop 1
    endif
    inv = inv/det
  end function inv3x3

  subroutine modulo1(x1)
    real(double) :: x1
    integer :: xint
    real(double) :: pluseps=1.0d-9

    if(x1 < 0.0d0) then

      xint = int(x1)
      x1 = x1 - xint + 1d0
      if(x1 == 1.0d0) x1 = 0d0
    else if(x1 >= 1.0d0) then

      xint = int(x1)
      x1 = x1 - xint
    endif

    if(x1 >= one .or. x1 < zero) then
      write(*,*) 'x1 is ',x1
      write(*,'(A)') 'err: must stop, x1 must be in [0,1) '
      stop 1
    endif

    if(abs(x1 - one) < pluseps) then

      x1 = zero
    endif
  end subroutine modulo1

  subroutine foldtofirstzone3(x)
    real(double) :: x(3)
    integer :: i
    do i = 1, 3
      call modulo1(x(i))
    enddo
  end subroutine foldtofirstzone3

  subroutine supercell_2_xyz(xyzu,filename,s)
    integer :: xyzu
    character(len=*) :: filename
    type(supercell) :: s
    integer :: i,n
    real(double) :: f(3),p(3)

    n = s%n
    open(unit=xyzu,file=trim(filename),status='replace')
    write(xyzu,'(I5)') n
    write(xyzu,'(A)') ' generated by VMD'
    do i = 1, n
      f(1:3) = s%at(i)%f(1:3)
      call foldtofirstzone3(f(1))
      call frac2abs(s%a,f(1),p(1))
      write(xyzu,'(A5,3F18.12,I4)') CatS(s%at(i)%z),p(1:3)
    enddo
    close(xyzu)

  end subroutine supercell_2_xyz

  subroutine supercell_2_siesta_fdf(siestau,filename,atmtype,zord,s)
    integer :: siestau
    character(len=*) :: filename
    type(supercell) :: s
    integer :: totat,i,ind,sizei,j,n,atmtype,zord(atmtype),groupsize(200)
    type(oneatom) :: outcell(maxnum)
    integer :: outcellatmtype(maxnum)
    integer,allocatable :: fmap(:)
    real(double) :: x(3),y(3)
    integer :: SortSwitch
    integer,parameter :: NoSorting=1,WithSorting=2
    integer,allocatable :: localbmap(:)

    n = s%n

    allocate(localbmap(n))

    allocate(fmap(n))
    ind = 0
    totat = 0

    do i = 1, atmtype
      sizei = 0
      do j = 1, n
        if(s%at(j)%z == zord(i)) then
          sizei = sizei + 1
          ind = ind + 1
          fmap(ind) = j
          outcell(ind)%z = zord(i)
          outcellatmtype(ind) = i
          if(trim(s%fracabs) == 'frac') then
            outcell(ind)%p(1:3) = s%at(j)%f(1:3)
          else if(trim(s%fracabs) == 'abs') then
            x = s%at(j)%f(1:3)
            call frac2abs(s%a,x,y)
            outcell(ind)%p(1:3) = y
          else
            write(*,'(A)') 'ocf is '//trim(s%fracabs)
            write(*,'(A)') 'err: ocf format not supported.'
            stop 1
          endif
        endif
      enddo
      groupsize(i) = sizei
      totat = totat + groupsize(i)
    enddo
    if(totat /= n) then
      write(*,*) 'totat,n=',totat,n
      write(*,'(A)') 'err: totat and n problem in supercell_2_siesta_fdf'
      stop 1
    endif
    if(ind /= n) then
      write(*,'(A,2I5)') 'ind,n = ',ind,n
      write(*,'(A)') 'err: new ind and n are not the same.'
      stop 1
    endif

    do i = 1, atmtype
      if(groupsize(i) < 1) then
        write(*,'(A,I4,A,I4)') 'i, groupsize(i)=',i,',',groupsize(i)
        write(*,'(A)') 'err: supercell_2_siesta_fdf: groupsize(i) is less than 1.'
        stop 1
      endif
    enddo

    do i = 1, n
      j = fmap(i)
      localbmap(j) = i
    enddo
    open(unit=siestau,file=trim(filename),status='replace')

    write(siestau,'(A,I5)') 'NumberOfAtoms ',n
    write(siestau,'(A,I5)') 'NumberOfSpecies ',atmtype
    write(siestau,'(A)') '%block ChemicalSpeciesLabel'
    do i = 1, atmtype
      write(siestau,'(I5,I5,A6)') i,zord(i),Cats(zord(i))
    enddo
    write(siestau,'(A)') '%endblock ChemicalSpeciesLabel'

    write(siestau,'(A)') 'LatticeConstant 1.0 Ang'
    write(siestau,'(A)') '%block LatticeVectors'
    write(siestau,'(3F18.12)') s%a(1:3,1)
    write(siestau,'(3F18.12)') s%a(1:3,2)
    write(siestau,'(3F18.12)') s%a(1:3,3)
    write(siestau,'(A)') '%endblock LatticeVectors'
    if(trim(s%fracabs) == 'frac') then
      write(siestau,'(A)') 'AtomicCoordinatesFormat  Fractional'
    else if(trim(s%fracabs) == 'abs') then
      write(siestau,'(A)') 'AtomicCoordinatesFormat  NotScaledCartesianAng'
    endif
    write(siestau,'(A)') '%block AtomicCoordinatesAndAtomicSpecies'

    SortSwitch=NoSorting

    do i = 1, s%n
      if(SortSwitch == WithSorting) then
        write(siestau,'(3F18.12,I4)') outcell(i)%p(1:3),outcellatmtype(i)
      else if(SortSwitch == NoSorting) then
        write(siestau,'(3F18.12,I4)') outcell(localbmap(i))%p(1:3),outcellatmtype(localbmap(i))
      else
        write(*,'(A)') 'err: logic error in supercell_2_siesta_fdf'
        stop 1
      endif
    enddo
    write(siestau,'(A)') '%endblock AtomicCoordinatesAndAtomicSpecies'
    close(siestau)
    deallocate(fmap)
    deallocate(localbmap)
  end subroutine supercell_2_siesta_fdf

  subroutine supercell_2_gulp(filename,atmtype,zord,s)
    character(len=*) :: filename
    type(supercell) :: s
    integer :: totat,i,ind,sizei,j,n,atmtype,zord(atmtype),groupsize(200)
    type(oneatom) :: outcell(maxnum)
    integer,parameter :: gulpu=50
    integer,allocatable :: fmap(:)
    character(len=200) :: commentline
    integer,allocatable :: localbmap(:)

    n = s%n
    allocate(fmap(n))
    allocate(localbmap(n))

    ind = 0
    totat = 0
    do i = 1, atmtype
      sizei = 0
      do j = 1, n
        if(s%at(j)%z == zord(i)) then
          sizei = sizei + 1
          ind = ind + 1
          fmap(ind) = j
          outcell(ind)%z = zord(i)
          outcell(ind)%p(1:3) = s%at(j)%f(1:3)
        endif
      enddo
      groupsize(i) = sizei
      totat = totat + groupsize(i)
    enddo
    if(totat /= n) then
      write(*,*) 'totat,n=',totat,n
      write(*,'(A)') 'err: totat and n problem in supercell_2_gulp'
      stop 1
    endif
    if(ind /= n) then
      write(*,'(A,2I5)') 'ind,n = ',ind,n
      write(*,'(A)') 'err: new ind and n are not the same.'
      stop 1
    endif

    do i = 1, atmtype
      if(groupsize(i) < 1) then
        write(*,'(A,I4,A,I4)') 'i, groupsize(i)=',i,',',groupsize(i)
        write(*,'(A)') 'err: supercell_2_gulp: groupsize(i) is less than 1.'
        stop 1
      endif
    enddo
    do i = 1, n
      j = fmap(i)
      localbmap(j) = i
    enddo
    open(unit=gulpu,file=trim(filename),status='replace')

    commentline='vectors'
    write(gulpu,'(A)') trim(commentline)
    write(gulpu,'(3F18.12)') s%a(1:3,1)
    write(gulpu,'(3F18.12)') s%a(1:3,2)
    write(gulpu,'(3F18.12)') s%a(1:3,3)
    commentline='fractional'
    write(gulpu,'(A)') trim(commentline)
    ind = 0
    do i = 1, atmtype
      do j = 1, groupsize(i)
        ind = ind + 1
       write(gulpu,'(A,A,3F18.12,A)') cats(zord(i)), ' core ',outcell(ind)%p(1:3)
      enddo
    enddo
    close(gulpu)
    deallocate(fmap)
    deallocate(localbmap)
  end subroutine supercell_2_gulp

  subroutine supercell_2_vasp(unitn,filename,atmtype,zord,s,writeforce,ATMSswitch)
    integer,optional :: writeforce
    integer,optional :: ATMSswitch
    integer :: unitn
    character(len=*) :: filename
    type(supercell) :: s
    integer :: atmtype,zord(atmtype)
    integer :: totat,i,ind,sizei,j,n
    integer :: groupsize(200)
    type(oneatom) :: outcell(maxnum)
    integer,allocatable :: fmap(:)
    character(len=DSL) :: commentline
    real(double),parameter :: maxlen=1000d0
    real(double) :: len6(6),frac(3),cartesian(3),num3(3)
    integer,allocatable :: localbmap(:)
    integer :: k

    n = s%n
    allocate(fmap(n))
    allocate(localbmap(n))

    ind = 0
    totat = 0
    do i = 1, atmtype
      sizei = 0
      do j = 1, n
        if(s%at(j)%z == zord(i)) then
          sizei = sizei + 1
          ind = ind + 1
          fmap(ind) = j
          outcell(ind)%z = zord(i)
          outcell(ind)%p(1:3) = s%at(j)%f(1:3)
          outcell(ind)%force(1:3) = s%at(j)%force(1:3)
        endif
      enddo
      groupsize(i) = sizei
      totat = totat + groupsize(i)
    enddo
    if(totat /= n) then
      write(*,*) 'totat,n=',totat,n
      write(*,'(A)') 'err: totat and n problem in supercell_2_vasp'
      stop 1
    endif
    if(ind /= n) then
      write(*,'(A,2I5)') 'ind,n = ',ind,n
      write(*,'(A)') 'err: new ind and n are not the same.'
      stop 1
    endif

    do i = 1, atmtype
      if(groupsize(i) < 1) then
        write(*,'(A,I4,A,I4)') 'i, groupsize(i)=',i,',',groupsize(i)
        write(*,'(A)') 'err: supercell_2_vasp: groupsize(i) is less than 1.'
        stop 1
      endif
    enddo
    do i = 1, n
      j = fmap(i)
      localbmap(j) = i
    enddo
    open(unit=unitn,file=trim(filename),status='replace')

    if(present(ATMSswitch)) then
      if(ATMSswitch == 1) then
        commentline = 'ATMS: '//trim(N2str(n))//' '//trim(N2str(atmtype))
        do i = 1, atmtype

          commentline = trim(commentline)//' '//trim(N2str(groupsize(i)))//' '//trim(Cats(zord(i)))
        enddo
        commentline = trim(commentline)//' SG '//trim(s%sg_label)
      endif
    else
      commentline = ''

      if(trim(s%sg_label) /= '' .and. trim(s%sg_label) /= 'undefined') then
        commentline = 'SG '//trim(s%sg_label)
      endif
    endif

    write(unitn,'(A)') trim(commentline)

    write(unitn,'(F6.3)') 1.0d0
    write(unitn,'(3F18.12)') s%a(1:3,1)
    write(unitn,'(3F18.12)') s%a(1:3,2)
    write(unitn,'(3F18.12)') s%a(1:3,3)

    commentline=''
    do i = 1, atmtype
      commentline = trim(commentline)//' '//trim(CATS(zord(i)))
    enddo
    write(unitn,'(A)') trim(commentline)
    write(unitn,'(10I4)') (groupsize(i),i=1,atmtype)
    if(trim(s%fracabs) == 'frac') then
      write(unitn,'(A)') 'Direct'
    else if(trim(s%fracabs) == 'abs') then
      write(unitn,'(A)') 'Cartesian'
    else
      write(*,'(A)') 'err: Must be either frac or abs.'
      stop 1
    endif

    ind = 0
    do i = 1, atmtype
      do j = 1, groupsize(i)
        ind = ind + 1
        frac(1:3) = outcell(ind)%p(1:3)

        call frac2abs(s%a,frac(1),cartesian(1))

        if(trim(s%fracabs) == 'frac') then
          num3(1:3) = frac(1:3)

          do k = 1, 3
            call modulo1(num3(k))
          enddo
        else if(trim(s%fracabs) == 'abs') then
          num3(1:3) = cartesian(1:3)
        endif
        if(.not. present(writeforce)) then
          write(unitn,'(3F18.12)') num3(1:3)
        else
          write(unitn,'(3F18.12,A,3F18.12)') num3(1:3),' | ',outcell(ind)%force(1:3)
        endif
      enddo
    enddo
    close(unitn)

    call getlenang(s%a,len6(1))
    deallocate(fmap)
    deallocate(localbmap)
  end subroutine supercell_2_vasp

  subroutine supercell_2_pwposcar(unitn,filename,atmtype,zord,s)
    integer :: unitn
    character(len=*) :: filename
    type(supercell) :: s
    integer :: totat,i,ind,sizei,j,n,atmtype,zord(atmtype),groupsize(200)
    type(oneatom) :: outcell(maxnum)
    integer,allocatable :: fmap(:)
    integer,allocatable :: localbmap(:)
    integer :: z
    character(len=2) :: zname

    if(s%n < 1) then
      write(*,'(A)') 'err: s%n is smaller than 1.'
      stop 1
    endif
    allocate(localbmap(s%n))

    if(trim(s%fracabs) == 'abs') then
      write(*,'(A)') 'ocf = '//trim(s%fracabs)
      write(*,*) 'the default is still in frac'
      write(*,*) 'err: not supported/implemented yet.'
      write(*,'(A)') 'err: check again.'
      stop 1
    endif
    n = s%n
    allocate(fmap(n))

    ind = 0
    totat = 0
    do i = 1, atmtype
      sizei = 0
      do j = 1, n
        if(s%at(j)%z == zord(i)) then
          sizei = sizei + 1
          ind = ind + 1
          fmap(ind) = j
          outcell(ind)%z = zord(i)
          outcell(ind)%p(1:3) = s%at(j)%f(1:3)
        endif
      enddo
      groupsize(i) = sizei
      totat = totat + groupsize(i)
    enddo
    if(totat /= n) then
      write(*,*) 'totat,n=',totat,n
      write(*,'(A)') 'err: totat and n problem in supercell_2_pwposcar'
      stop 1
    endif
    if(ind /= n) then
      write(*,'(A,2I5)') 'ind,n = ',ind,n
      write(*,'(A)') 'err: new ind and n are not the same.'
      stop 1
    endif

    do i = 1, atmtype
      if(groupsize(i) < 1) then
        write(*,'(A,I4,A,I4)') 'i, groupsize(i)=',i,',',groupsize(i)
        write(*,'(A)') 'err: supercell_2_pwposcar: groupsize(i) is less than 1.'
        stop 1
      endif
    enddo
    do i = 1, n
      j = fmap(i)
      localbmap(j) = i
    enddo

    open(unit=unitn,file=trim(filename)//'1',status='replace')
    write(unitn,'(A)') 'ibrav=0'
    write(unitn,'(A)') 'nat='//trim(N2str(n))
    write(unitn,'(A)') 'ntyp='//trim(N2str(atmtype))
    write(unitn,'(A)') 'celldm(1)=1.0'
    close(unit=unitn)

    open(unit=unitn,file=trim(filename)//'2',status='replace')
    write(unitn,'(A)') 'ATOMIC_SPECIES'
    do i = 1, atmtype
      z = zord(i)
      write(unitn,'(A,A,F18.12,A,A)') trim(CatS(z)),'   ',MASSOFATOM(z),'   ',trim(CatS(z))//'.pp'
    enddo
    write(unitn,'(A)') 'CELL_PARAMETERS'
    write(unitn,'(3F18.12)') s%a(1:3,1)
    write(unitn,'(3F18.12)') s%a(1:3,2)
    write(unitn,'(3F18.12)') s%a(1:3,3)
    write(unitn,'(A)') 'ATOMIC_POSITIONS {crystal}'
    ind = 0
    do i = 1, atmtype
      do j = 1, groupsize(i)
        ind = ind + 1
        zname = CatS(outcell(ind)%z)
        zname = adjustl(zname)
        write(unitn,'(A2,A,3F18.12)') zname,'   ',outcell(ind)%p(1:3)
      enddo
    enddo
    close(unitn)
    deallocate(fmap)
    deallocate(localbmap)
  end subroutine supercell_2_pwposcar

  subroutine supercell_2_fleur(unitn,filename,atmtype,zord,s)
    integer :: unitn
    character(len=*) :: filename
    type(supercell) :: s
    integer :: totat,i,ind,sizei,j,n,atmtype,zord(atmtype),groupsize(200)
    type(oneatom) :: outcell(maxnum)
    integer,allocatable :: fmap(:)
    integer,allocatable :: localbmap(:)

    if(trim(s%fracabs) == 'abs') then
      write(*,'(A)') 'ocf = '//trim(s%fracabs)
      write(*,*) 'the default is still in frac'
      write(*,*) 'not supported/implemented yet.'
      write(*,'(A)') 'err: check again.'
      stop 1
    endif
    n = s%n
    allocate(fmap(n))
    allocate(localbmap(n))

    ind = 0
    totat = 0
    do i = 1, atmtype
      sizei = 0
      do j = 1, n
        if(s%at(j)%z == zord(i)) then
          sizei = sizei + 1
          ind = ind + 1
          fmap(ind) = j
          outcell(ind)%z = zord(i)
          outcell(ind)%p(1:3) = s%at(j)%f(1:3)
        endif
      enddo
      groupsize(i) = sizei
      totat = totat + groupsize(i)
    enddo
    if(totat /= n) then
      write(*,*) 'totat,n=',totat,n
      write(*,'(A)') 'err: totat and n problem in supercell_2_pwposcar'
      stop 1
    endif
    if(ind /= n) then
      write(*,'(A,2I5)') 'ind,n = ',ind,n
      write(*,'(A)') 'err: new ind and n are not the same.'
      stop 1
    endif

    do i = 1, atmtype
      if(groupsize(i) < 1) then
        write(*,'(A,I4,A,I4)') 'i, groupsize(i)=',i,',',groupsize(i)
        write(*,'(A)') 'err: supercell_2_pwposcar: groupsize(i) is less than 1.'
        stop 1
      endif
    enddo
    do i = 1, n
      j = fmap(i)
      localbmap(j) = i
    enddo

    open(unit=unitn,file=trim(filename),status='replace')
    write(unitn,'(A)') 'CELL_PARAMETERS'
    write(unitn,'(3F18.12)') s%a(1:3,1)
    write(unitn,'(3F18.12)') s%a(1:3,2)
    write(unitn,'(3F18.12)') s%a(1:3,3)
    write(unitn,*)
    write(unitn,'(A)') trim(N2str(s%n))
    ind = 0
    do i = 1, atmtype
      do j = 1, groupsize(i)
        ind = ind + 1
        write(unitn,'(I4,3F18.12)') s%at(ind)%z,outcell(ind)%p(1:3)
      enddo
    enddo
    close(unitn)
    deallocate(fmap)
    deallocate(localbmap)
  end subroutine supercell_2_fleur

  subroutine supercell_2_xcrysden(xcrysdenfile,s)
    type(supercell) :: s
    character(len=*) :: xcrysdenfile
    real(double) :: x(3)
    integer,parameter :: xcrysdenu=50
    integer :: i

    open(unit=xcrysdenu,file=trim(xcrysdenfile),status='replace')
    write(xcrysdenu,'(A)') 'CRYSTAL'
    write(xcrysdenu,'(A)') 'PRIMVEC'
    do i = 1, 3
      write(xcrysdenu,'(3f20.10)') s%a(1:3,i)
    enddo
    write(xcrysdenu,'(A)') 'CONVVEC'
    do i = 1, 3
      write(xcrysdenu,'(3f20.10)') s%a(1:3,i)
    enddo
    write(xcrysdenu,'(A)') 'PRIMCOORD'
    if(s%n < 1) then
      write(*,*) 's%n = ',s%n
      write(*,'(A)') 'err: invalid n in supercell_2_xcrysden'
      stop 1
    endif
    write(xcrysdenu,'(2I4)') s%n, 1
    do i = 1, s%n
      call frac2abs(s%a,s%at(i)%f(1:3),x(1:3))
      write(xcrysdenu,'(I4,3F20.10)') s%at(i)%z,x(1:3)
    enddo
    close(xcrysdenu)
  end subroutine supercell_2_xcrysden

  subroutine onespace(ns,comment1,comment)
    integer :: ns
    character(len=ns) :: comment1,comment
    integer :: i,j

    comment = 's'

    i = 1
    j = 1

    loop2: do
      loop1: do
        if(comment1(i:i) /= ' ') exit loop1
        i = i + 1
        if(i == ns) exit loop2
      enddo loop1

      loop3: do
        if(comment1(i:i) == ' ') exit loop3
        comment(j:j) = comment1(i:i)
        j = j + 1
        i = i + 1
      enddo loop3
      comment(j:j) = ' '
      j = j + 1
    enddo loop2
  end subroutine onespace

  subroutine assign_str_fr_struct_out(inu,s)
    integer :: inu
    type(supercell) :: s
    integer :: i

    do i = 1, 3
      read(inu,*) s%a(1:3,i)
    enddo
    read(inu,*) s%n
    do i = 1, s%n
      read(inu,*) s%at(i)%gr, s%at(i)%z,s%at(i)%f(1:3)
    enddo
    call real2recip(s%a,s%b)
    close(inu)
  end subroutine assign_str_fr_struct_out

  subroutine assign_str_fr_xsf(controlu,sc)
    integer :: i,controlu
    type(supercell) :: sc
    character(LEN=DSL) :: tmp

    read(controlu,*) tmp
    if(trim(tmp) /= 'CRYSTAL') then
      write(*,'(A)') 'err: Keyword CRYSTAL not found.'
      stop 1
    endif
    read(controlu,*) tmp
    if(trim(tmp) /= 'PRIMVEC') then
      write(*,'(A)') 'err: Keyword PRIMVEC not found.'
      stop 1
    endif

    read(controlu,*) sc%a(1:3,1)
    read(controlu,*) sc%a(1:3,2)
    read(controlu,*) sc%a(1:3,3)
    read(controlu,*) tmp
    if(trim(tmp) /= 'CONVVEC') then
      write(*,'(A)') 'err: Keyword CONVVEC not found.'
      stop 1
    endif

    read(controlu,*) sc%a(1:3,1)
    read(controlu,*) sc%a(1:3,2)
    read(controlu,*) sc%a(1:3,3)
    call real2recip(sc%a,sc%b)
    read(controlu,*) tmp
    if(trim(tmp) /= 'PRIMCOORD') then
      write(*,'(A)') 'err: Keyword PRIMCOORD not found.'
      stop 1
    endif
    read(controlu,*) sc%n
    do i = 1, sc%n
      read(controlu,*) sc%at(i)%z, sc%at(i)%ac
      call abs2frac(sc%b,sc%at(i)%ac,sc%at(i)%f)
    enddo
    sc%sg_label = '1-P1'
  end subroutine assign_str_fr_xsf

  subroutine assign_str_fr_poscar(controlu,sc)
    integer :: controlu
    type(supercell) :: sc
    real(double) :: totalscaling,vol
    integer :: i,j,s,nsp,v,z,ind
    character(LEN=DSL) :: comment1,comment2,tmp1,tmp2,valuestr,targetlabel
    real(double) :: pos(3)
    integer,allocatable :: numarr(:)
    integer :: totn,k,labellen,SGnumber,loc
    logical :: done,found,puredigits
    integer :: nvalid
    real(double) :: density

    read(controlu,DSF) comment1

    call onespace(DSL,comment1,comment2)

    sc%commentline = trim(comment2)
    write(*,'(A,A,A)') 'POSCAR comment is [',trim(comment2),'].'
    s = scan(comment2,' ')

    if(comment2(1:s-1) == 'ATMS:') then
      tmp1= comment2(s+1:DSL)

      read(tmp1,*) v
      sc%n = v
      if(sc%n < 1) then
        write(*,*) 'sc%n = ',sc%n
        write(*,'(A)') 'err: check sc%n in assign_str_fr_poscar'
        stop 1
      endif
      s = scan(tmp1,' ')
      tmp1 = tmp1(s+1:DSL)

      read(tmp1,*) v
      sc%nsp = v
      if(sc%nsp < 1) then
        write(*,*) 'sc%nsp = ',sc%nsp
        write(*,'(A)') 'err: check sc%nsp in assign_str_fr_poscar.'
        stop 1
      elseif(sc%nsp > zordmax) then
        write(*,*) 'sc%nsp,zordmax=',sc%nsp,zordmax
        write(*,'(A)') 'err: check sc%nsp, i.e., the number of species/elements'
        stop 1
      endif
      s = scan(tmp1,' ')
      tmp1 = tmp1(s+1:DSL)
      nsp = v
      if(sc%nsp > zordmax) then
        write(*,*) 'sc%nsp = ',sc%nsp
        write(*,*) 'zordmax = ',zordmax
        write(*,'(A)') 'err: out of bound.'
        stop 1
      endif

      do i = 1, nsp

        read(tmp1,*) v

        sc%n_per_species(i) = v

        s = scan(tmp1,' ')
        tmp1 = tmp1(s+1:DSL)

        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        found = .false.
        do j = 1, zordmax
          if(found) exit
          if(trim(valuestr) == trim(Ats(j)) .or. trim(valuestr) == trim(CAts(j)) ) then
            z = j
            found = .true.
          endif
        enddo
        if(.not. found) then
          write(*,*) 'valuestr = ',trim(valuestr)
          write(*,'(A)') 'err: this valuestr is not found in the periodic table.'
          stop 1
        endif
        sc%zarr(i) = z
        tmp1 = tmp1(s+1:DSL)
      enddo

      ind = 0
      do j = 1, sc%nsp
        ind = ind + sc%n_per_species(j)
      enddo

      if(ind /= sc%n) then
        write(*,*) 'ind,sc%n=',ind,sc%n
        write(*,'(A)') 'err: missing atoms in the assign_str_fr_poscar.'
        stop 1
      endif

      sc%sg_label = ''
      if(trim(tmp1) /= '') then

        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        if(trim(valuestr) /= 'SG') then
          write(*,*) 'valuestr is '//trim(valuestr)
          write(*,*) 'The first optional parameter must be SG.'
          write(*,'(A)') 'err: Check poscar.'
          stop 1
        endif
        tmp1 = tmp1(s+1:DSL)
        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        sc%sg_label = trim(valuestr)
      else
        write(*,*) 'No SG label'
      endif

      read(controlu,*) totalscaling
      read(controlu,*) sc%a(1:3,1)
      read(controlu,*) sc%a(1:3,2)
      read(controlu,*) sc%a(1:3,3)

      allocate(numarr(sc%nsp))

      read(controlu,DSF) tmp1
      write(*,*) 'tmp1 is .'//trim(tmp1)//'.'
      call onespace(DSL,tmp1,tmp2)
      s = scan(tmp2,'0123456789')

      if(s == 0) then

        write(*,'(A)') 'We are going to process the species line |'//trim(tmp2)//'|'

        nsp = num_of_strings(tmp2)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif
        do i = 1, sc%nsp
          s = scan(tmp2,' ')
          valuestr = tmp2(1:s-1)
          write(*,*) 'valuestr is '//trim(valuestr)
          found = .false.
          do j = 1, zordmax
            if(found) exit
            if(trim(valuestr) == trim(CATS(j))) then
              z = j
              if(z /= sc%zarr(i)) then
                write(*,*) 'z, sc%zarr(i) = ',z,sc%zarr(i)
                write(*,'(A)') 'err: Inconsistency.'
                stop 1
              endif
              found = .true.
            endif
          enddo
          if( .not. found) then
            write(*,'(A)') 'err: not found.'
            stop 1
          endif
          tmp2 = tmp2(s+1:DSL)
        enddo

        read(controlu,DSF) tmp2
        nsp = num_of_integers(tmp2)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif

        read(tmp2,*) numarr(1:sc%nsp)
      else

        nsp = num_of_integers(tmp1)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif
        read(tmp1,*) numarr(1:sc%nsp)
      endif

      do i = 1, sc%nsp
        if(numarr(i) /= sc%n_per_species(i)) then
          write(*,*) 'species: i=',i
          write(*,*) 'numarr(i) = '//trim(N2str(numarr(i)))//', sc%n_per_species(i) = '//trim(N2str(sc%n_per_species(i)))
          write(*,'(A)') 'err: numarr(i) and sc%n_per_species(i) are not consistent'
          stop 1
        endif
      enddo
      deallocate(numarr)

    else
      if(comment2(1:s-1) == 'SG') then

        tmp2 = comment2(s+1:DSL)
        targetlabel = trim(tmp2)

        write(*,*) 'targetlabel = '//trim(targetlabel)

        puredigits = .true.
        labellen = len(trim(targetlabel))
        do i = 1, labellen
          if(.not. puredigits) then
            exit
          endif
         loc = scan(targetlabel(i:i),'0123456789')
         if(loc /= 1) then
           puredigits = .false.
         endif
        enddo

        if(puredigits)  then
          SGnumber = str2N(trim(targetlabel))
          write(*,*) 'SGnumber is ',SGnumber
          if(SGnumber >= 1 .and. SGnumber <= 230) then
            targetlabel = trim(SGbase(1,SGnumber)%sglabel)
            write(*,*) 'Fully expand the default SG label to the full SG label:'
            write(*,*) 'SGnumber= ',SGnumber
            write(*,*) 'Full targetlabel= ',trim(targetlabel)
          else
            write(*,*) 'err: must be between 1 and 230 inclusive. impossible. '
            stop 1
          endif
        endif
        sc%sg_label = trim(targetlabel)
        write(*,'(A)') 'Possibly adjusted sg_label is '//trim(sc%sg_label)
      else
        write(*,'(A)') 'WARNING WARNING: No ATMS:, no SG, hence we assume SG 1'
        sc%sg_label = '1-P1'
      endif

      read(controlu,*) totalscaling
      read(controlu,*) sc%a(1:3,1)
      read(controlu,*) sc%a(1:3,2)
      read(controlu,*) sc%a(1:3,3)

      read(controlu,DSF) tmp1
      call onespace(DSL,tmp1,tmp2)

      s = scan(tmp2,'0123456789')

      if(s == 0) then

        nsp = 0
        done = .false.
        do
          if(done) exit
          s = scan(tmp2,' ')
          valuestr = tmp2(1:s-1)

          found = .false.
          do j = 1, zordmax
            if(found) exit
            if(trim(valuestr) == trim(CATS(j))) then
              nsp = nsp + 1
              sc%zarr(nsp) = j
              found = .true.
            endif
          enddo
          if( .not. found) then
            write(*,'(A)') 'err: not found.'
            stop 1
          endif
          tmp2 = tmp2(s+1:DSL)
          if(len(trim(tmp2)) == 0) then
            write(*,'(A)') 'nsp = '//trim(N2str(nsp))
            sc%nsp = nsp
            done = .true.
          endif
        enddo

        read(controlu,DSF) tmp1

        nvalid = num_of_integers(tmp1)
        if(nsp /= nvalid) then
          write(*,*) 'nsp,nvalid=',nsp,nvalid
          write(*,*) 'err: inconsistency in species line and species line.'
          stop 1
        endif
        read(tmp1,*) sc%n_per_species(1:nsp)

      else
        write(*,'(A)') 'err: Since this is without ATMs, we must have species line.'
        stop 1
      endif
    endif

    totn = 0
    do j = 1, sc%nsp
      v = sc%n_per_species(j)
      do k = 1, v
        totn = totn + 1
        sc%at(totn)%z = sc%zarr(j)
        sc%at(totn)%mass = massofa(sc%zarr(j))
      enddo
    enddo
    write(*,'(A)') 'totn = '//trim(N2str(totn))
    sc%n = totn

    sc%a = sc%a*totalscaling

    call real2recip(sc%a,sc%b)

    write(*,'(A,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10)') &
      'Real latt:', sc%a(1:3,1),'    | ',vecmag3(sc%a(1:3,1)),sc%a(1:3,2),'    | ',vecmag3(sc%a(1:3,2)),sc%a(1:3,3),'    | ',vecmag3(sc%a(1:3,3))
    vol = det3(sc%a(1,1))
    if(vol .le. 0) then
      write(*,'(A)') 'err: vol is not positive.'
      stop 1
    endif

    sc%vol = vol

    write(*,'(A,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10)') &
      'Reciprocal latt:', sc%b(1:3,1),'    | ',vecmag3(sc%b(1:3,1)),sc%b(1:3,2),'    | ',vecmag3(sc%b(1:3,2)),sc%b(1:3,3),'    | ',vecmag3(sc%b(1:3,3))

    read(controlu,*) tmp1
    if(trim(tmp1) == 'Selective') then
      read(controlu,*) tmp1
    endif
    if(trim(tmp1) /= 'Direct' .and. trim(tmp1) /= 'direct' .and. trim(tmp1) /= 'Cartesian' ) then
      write(*,'(A)') 'err: dummy should be Direct, direct, or Cartesian'
      stop 1
    endif

    do i = 1, sc%n
      sc%at(i)%force(1:3) = (/zero,zero,zero/)
    enddo
    do i = 1, sc%n
      if(trim(tmp1) == 'Direct' .or. trim(tmp1) == 'direct') then
        read(controlu,DSF) comment1
        call onespace(DSL,comment1,comment2)
        read(comment2,*) sc%at(i)%f(1:3)

        call frac2abs(sc%a(1,1),sc%at(i)%f(1),sc%at(i)%ac(1))
        s = scan(comment2,'|')
        if(s > 0) then
          tmp2= comment2(s+1:DSL)

          read(tmp2,*) sc%at(i)%force(1:3)
          write(*,'(A,3F20.10)') 'sc%at(i)%force(1:3) = ',sc%at(i)%force(1:3)
        endif
      else if(trim(tmp1) == 'Cartesian') then

        read(controlu,DSF) comment1
        call onespace(DSL,comment1,comment2)
        read(comment2,*) pos(1:3)

        call abs2frac(sc%b(1,1),pos(1),sc%at(i)%f(1))

        call frac2abs(sc%a(1,1),sc%at(i)%f(1),sc%at(i)%ac(1))
        s = scan(comment2,'|')
        if(s > 0) then
          tmp2= comment2(s+1:DSL)

          read(tmp2,*) sc%at(i)%force(1:3)
          write(*,'(A,3F20.10)') 'sc%at(i)%force(1:3) = ',sc%at(i)%force(1:3)
        endif
      endif
    enddo
    call getlenang(sc%a,sc%la(1))
    write(*,'(A)') 'la(1:6)='
    do i = 1, 6
      write(*,'(3F15.10)') sc%la(i)
    enddo
    write(*,*) 'Volume of the supercell= ',vol
    density = crys_density(sc)
    sc%density = density
    write(*,*) 'Density is ',density,' g/cm^3'

  end subroutine assign_str_fr_poscar

  subroutine assign_str_fr_fdf(inu,s)
    integer :: inu,natom,nspecies,sp,i,j
    type(supercell) :: s
    character(len=DSL) :: str1,str2,str3
    real(double) :: x(3)
    character(len=DSL) :: icf

    read(inu,*) str1,str2
    if(trim(str1) /= 'NumberOfAtoms') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be NumberofAtoms'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    else
      read(str2,'(I8)') natom
      write(*,*) 'natom = ',natom
      s%n = natom
    endif
    read(inu,*) str1,str2
    if(trim(str1) /= 'NumberOfSpecies') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be NumberofSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    else
      read(str2,'(I8)') nspecies
      write(*,*) 'nspecies = ',nspecies
      s%nsp = nspecies
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'ChemicalSpeciesLabel') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be ChemicalSpeciesLabel'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    do i = 1, nspecies
      read(inu,*) j,sp
      if(j/=i) then
        write(*,*) 'j,i=',j,i
        write(*,*) 'sp = ',sp
        write(*,'(A)') 'err: something wrong with the species numbering ?'
        stop 1
      else
        s%zarr(i) = sp
        write(*,*) 's%zarr(',i,')=',sp
      endif
    enddo
    read(inu,*) str1,str2
    if(trim(str2) /= 'ChemicalSpeciesLabel') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be ChemicalSpeciesLabel'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2,str3
    if(trim(str1) /= 'LatticeConstant') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be LatticeConstant'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    if(trim(str2) /= '1.0') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be 1.0'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    if(trim(str3) /= 'Ang') then
      write(*,*) 'str3 = ',trim(str3)
      write(*,*) 'but str3 must be Ang'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'LatticeVectors') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be  LatticeVectors'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) s%a(1:3,1)
    read(inu,*) s%a(1:3,2)
    read(inu,*) s%a(1:3,3)
    write(*,*) 's%a is '
    write(*,*) s%a(1:3,1)
    write(*,*) s%a(1:3,2)
    write(*,*) s%a(1:3,3)
    call real2recip(s%a,s%b)
    read(inu,*) str1,str2
    if(trim(str2) /= 'LatticeVectors') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be LatticeVectors'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) == 'Fractional') then
      icf = 'frac'
    else if(trim(str2) == 'NotScaledCartesianAng') then
      icf = 'abs'
    else
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be either Fractional or NotScaledCartesianAng'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'AtomicCoordinatesAndAtomicSpecies') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be AtomicCoordinatesAndAtomicSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    do i = 1, natom
      read(inu,*) x(1:3),sp
      if(trim(icf) == 'frac') then
        s%at(i)%f(1:3) = x(1:3)
      else if(trim(icf) == 'abs') then
        call abs2frac(s%b,x,s%at(i)%f)
      else
        write(*,'(A)') 'err: frac and abs problem.'
        stop 1
      endif
      s%at(i)%z = s%zarr(sp)

    enddo
    read(inu,*) str1,str2
    if(trim(str2) /= 'AtomicCoordinatesAndAtomicSpecies') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be AtomicCoordinatesAndAtomicSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
  end subroutine assign_str_fr_fdf

  subroutine supercell_2_cif(filename,atmtype,zord,s)
    character(len=*) :: filename
    type(supercell) :: s
    integer :: z,i,ind,sizei,j,n,atmtype,zord(zordmax),&
      groupsize(200)
    real(double) :: p(3)
    type(oneatom) :: outcell(maxnum)
    integer,parameter :: cifu=50
    real(double) :: L(6)
    integer,allocatable :: fmap(:)
    integer,allocatable :: localbmap(:)

    n = s%n
    allocate(fmap(n))
    allocate(localbmap(n))

    ind = 0
    do i = 1, atmtype
      sizei = 0
      do j = 1, n
        if(s%at(j)%z == zord(i)) then
          sizei = sizei + 1
          ind = ind + 1
          fmap(ind) = j
          outcell(ind)%z = zord(i)
          outcell(ind)%p(1:3) = s%at(j)%f(1:3)
        endif
      enddo
      groupsize(i) = sizei
    enddo
    do i = 1, n
      j = fmap(i)
      localbmap(j) = i
    enddo
    if(ind /= n) then
      write(*,'(A,2I5)') 'ind,n = ',ind,n
      write(*,'(A)') 'err: new ind and n are not the same.'
      stop 1
    endif
    call getlenang(s%a,L(1))
    open(unit=cifu,file=trim(filename),status='replace')
    write(cifu,'(A)') 'data_genericname'
    write(cifu,'(A)') '_audit_creation_date              2005-10-31'

    write(cifu,'(A)') '_symmetry_space_group_name_H-M    ''P1'''
    write(cifu,'(A)') '_symmetry_Int_Tables_number       1'
    write(cifu,'(A)') '_symmetry_cell_setting            triclinic'
    write(cifu,'(A)') 'loop_'
    write(cifu,'(A)') '_symmetry_equiv_pos_as_xyz'
    write(cifu,'(A)') '  x,y,z'
    write(cifu,'(A,F12.7)') '_cell_length_a',L(1)
    write(cifu,'(A,F12.7)') '_cell_length_b',L(2)
    write(cifu,'(A,F12.7)') '_cell_length_c',L(3)
    write(cifu,'(A,F12.7)') '_cell_angle_alpha',L(4)
    write(cifu,'(A,F12.7)') '_cell_angle_beta',L(5)
    write(cifu,'(A,F12.7)') '_cell_angle_gamma',L(6)
    write(cifu,'(A)') 'loop_'
    write(cifu,'(A)') '_atom_site_label'
    write(cifu,'(A)') '_atom_site_type_symbol'
    write(cifu,'(A)') '_atom_site_fract_x'
    write(cifu,'(A)') '_atom_site_fract_y'
    write(cifu,'(A)') '_atom_site_fract_z'
    write(cifu,'(A)') '_atom_site_occupancy'

    ind = 0
    do i = 1, atmtype
      do j = 1, groupsize(i)
        ind = ind + 1
        z = outcell(ind)%z
        p = outcell(ind)%p(1:3)
        write(cifu,'(A,A4,3F12.7,A,F4.1)') &
          trim(CAts(z))//trim(N2str(ind)),CAts(z),p(1:3),'   ',1.0
      enddo
    enddo
    close(cifu)
    deallocate(fmap)
    deallocate(localbmap)
  end subroutine supercell_2_cif

  function massofa(z)
    integer :: z
    real(double) :: massofa
    if(z < 0 .or. z > zordmax) then
      write(*,*) 'z = ',z
      write(*,'(A)') 'massofa not defined yet.'
      write(*,*) 'err: massofa problem'
      stop 1
    endif
    massofa = MASSOFATOM(z)
  end function massofa

  function crys_density(s)
    type(supercell) :: s
    real(double) :: mass,crys_density
    integer :: i

    mass = 0.0d0
    do i = 1, s%n
      mass = mass + massofatom(s%at(i)%z)
    enddo
    crys_density = (mass*AMU/(det3(s%a(1:3,1:3))*1.d-30))*(1.0d3/1.0d6)
  end function crys_density

  subroutine dbl_sort(n,v,ind,kflag)
    integer,parameter :: maxrecur=50
    integer :: i,n,segment,tind,ind(n),eindex(maxrecur),bindex(maxrecur),bi,ei,kflag,rightindex,leftindex
    real(double) :: v(n),vref,t
    logical :: foundright,foundleft,cross,inseq

    if(kflag /= 1 .and. kflag /= -1) then
      write(*,*) 'kflag in dbl_sort = ',kflag
      write(*,'(A)') 'err: kflag should be 1 (ascending) or -1 (descending).'
      stop 1
    endif

    if(n < 1) then
      write(*,*) 'dbl_sort, n = ',n
      write(*,'(A)') 'err: wrong array size for sorting.'
      stop 1
    endif
    if(n == 1) then
      return
    endif

    inseq = .true.
    i = 2
    do
      if(i > n) then
        exit
      endif
      if(.not. inseq) then
        exit
      endif
      if(kflag == 1) then
        inseq = v(i) >= v(i-1)
      else if(kflag == -1) then
        inseq = v(i) <= v(i-1)
      endif
      i = i + 1
    enddo
    if(inseq) then
     return
    endif

    if(kflag == -1) then
      do i = 1, n
        v(i) = -v(i)
      enddo
    endif

    segment = 1
    bindex(1) = 1
    eindex(1) = n
    do

      if(segment == 0) exit

      bi = bindex(segment)
      ei = eindex(segment)
      vref = v(bi)
      rightindex = ei
      leftindex = bi+1
      cross = .false.
      do
        if(cross) then
          exit
        endif
        foundright = .false.
        do
          if(foundright) exit
          if(v(rightindex) >= vref) then
            rightindex = rightindex - 1
            if(rightindex < bi) then
              foundright = .true.
            endif
          else
            foundright = .true.
          endif
        enddo
        foundleft = .false.
        do
          if(foundleft) exit
          if(v(leftindex) <= vref) then
            leftindex = leftindex + 1
            if(leftindex > ei) then
              foundleft = .true.
            endif
          else
            foundleft = .true.
          endif
        enddo
        if(leftindex > rightindex) then
          cross = .true.

          if(rightindex < bi) then

            segment = segment - 1

            if(ei-bi > 1) then
              segment = segment + 1
              bindex(segment) = bi+1
              eindex(segment) = ei
            endif
          else

            if(v(rightindex) < vref) then
              v(bi) = v(rightindex)
              v(rightindex) = vref
              tind = ind(bi)
              ind(bi) = ind(rightindex)
              ind(rightindex) = tind
            endif

            segment = segment - 1

            if(rightindex - bi > 1) then
              segment = segment + 1
              bindex(segment) = bi
              eindex(segment) = rightindex - 1
            endif
            if(ei-rightindex > 1) then
              segment = segment + 1
              if(segment > maxrecur) then
                write(*,*) 'maxrecur = ',maxrecur
                write(*,'(A)') 'maxrecur is not large enough'
              endif
              bindex(segment) = rightindex + 1
              eindex(segment) = ei
            endif
          endif
        else
          t = v(leftindex)
          v(leftindex) = v(rightindex)
          v(rightindex) = t
          tind = ind(leftindex)
          ind(leftindex) = ind(rightindex)
          ind(rightindex) = tind
          leftindex = leftindex+1
          rightindex = rightindex-1
        endif
      enddo
    enddo
    do i = 2, n
      if(v(i) < v(i-1)) then
        write(*,*) 'v(i),v(i-1)=',v(i),v(i-1)
        write(*,'(A)') 'err: v(i) must be greater or equal to v(i-1).'
        stop 1
      endif
    enddo
    if(kflag == -1) then
      do i = 1, n
        v(i) = -v(i)
      enddo
    endif
  end subroutine dbl_sort

  subroutine nearest_neighbor_tab(fu,uc,n1,n2,n3,nmax)
    integer :: fu,n1,n2,n3
    type(supercell) :: uc
    integer,parameter :: outputlen=1000
    character(len=outputlen) :: Line,tmpstr

    integer :: nmax
    real(double),allocatable :: distance(:)
    integer,allocatable :: atomindex(:),indarr(:)
    real(double) :: dis,nf(3),p1(3),p2(3)
    integer :: ind,i,j,k,m,n,z,arrlen
    real(double),parameter :: cut_dis=10.0d0

    if(n1 < 1 .or. n2  < 1 .or. n3 < 1) then
      write(*,*) 'n1, n2, n3 = ',n1,n2,n3
      write(*,'(A)') 'err: wrong input in nearest_neighbor.'
      stop 1
    endif

    allocate(distance(nmax))
    allocate(atomindex(nmax))
    allocate(indarr(nmax))

    write(fu,'(A)') 'Nearest neighbor table:'
    do i = 1, uc%n
      call frac2abs(uc%a,uc%at(i)%f,p1)

      do j = 1, nmax
        distance(j) = 1.0d100
      enddo

      ind = 0
      do j = -n1, n1
        do k = -n2, n2
          do m = -n3, n3
            do n = 1, uc%n
              nf = uc%at(n)%f + (/j,k,m/)
              call frac2abs(uc%a,nf,p2)
              dis = vecmag3(p1-p2)

              if(dis < cut_dis) then
                ind = ind + 1
                if(ind > nmax) then
                  write(*,*) 'nmax is not large enough to hold cut_dis'
                  write(*,*) 'nmax,cut_dis= ',nmax,cut_dis
                  write(*,'(A)') 'err: out of bound.'
                  stop 1
                endif
                distance(ind) = dis
                atomindex(ind) = n
              endif
            enddo
          enddo
        enddo
      enddo

      arrlen = ind
      write(*,*) 'total number of distances recorded = ',arrlen
      do j = 1, arrlen
        indarr(j) = j
      enddo
      call dbl_sort(arrlen,distance(1),indarr(1),1)

      z = i
      write(tmpstr,'(I0)') z
      Line = 'atm '//trim(tmpstr)
      write(tmpstr,'(A)') '('//trim(cats(uc%at(z)%z))//')'
      Line = trim(Line)//trim(tmpstr)//' |'

      do j = 2, arrlen
        k = indarr(j)
        write(tmpstr,'(I4)') atomindex(k)
        Line = trim(Line)//' '//trim(tmpstr)

        write(tmpstr,'(A)') trim(cats(uc%at(atomindex(k))%z))
        Line = trim(Line)//'('//trim(tmpstr)//')'
        write(tmpstr,'(F6.3)') distance(j)
        Line = trim(Line)//' '//trim(tmpstr)
        Line = trim(Line)//'|'
      enddo
      write(*,'(A)') trim(Line)
      write(fu,'(A)') trim(Line)
    enddo
    deallocate(distance)
    deallocate(atomindex)
    deallocate(indarr)
  end subroutine nearest_neighbor_tab

  subroutine get_nsp_zord(t,atmtype,zord)
    type(supercell) :: t
    integer :: atmtype,zord(zordmax)
    integer :: ind,i,Z
    logical :: found

    atmtype = 0
    do i = 1, t%n
      z = t%at(i)%z
      found = .false.
      ind = 0
      do
        if (found) exit
        if (ind > atmtype) exit
        ind = ind + 1
        if(z == zord(ind)) then
          found = .true.
        endif
      enddo

      if (.not. found) then
        atmtype = atmtype + 1
        zord(atmtype) = z
      endif
    enddo
    if(atmtype == 0) then
      write(*,'(A)') 'err: atmtype is zero.'
      stop 1
    endif
    write(*,*) 'atmtype = ',atmtype
    write(*,*) zord(1:atmtype)
  end subroutine get_nsp_zord

  subroutine cp_cell_with_scaling(oldstr,newstr,factor)
    type(supercell) :: oldstr,newstr
    real(double) :: factor
    integer :: i

    newstr%a(1:3,1) = oldstr%a(1:3,1)*factor
    newstr%a(1:3,2) = oldstr%a(1:3,2)*factor
    newstr%a(1:3,3) = oldstr%a(1:3,3)*factor
    newstr%n = oldstr%n
    newstr%nsp = oldstr%nsp
    do i = 1, oldstr%n
      newstr%at(i)%z = oldstr%at(i)%z
      newstr%at(i)%f = oldstr%at(i)%f
      call frac2abs(newstr%a,newstr%at(i)%f,newstr%at(i)%ac)

      newstr%at(i)%chg = oldstr%at(i)%chg
    enddo
  end subroutine cp_cell_with_scaling

  subroutine read_struc(inputu,inputfile,filestem,inputformat,s)
    integer :: inputu,length,ind
    character(len=*) :: inputfile,filestem,inputformat
    type(supercell) :: s

    length = len_trim(inputfile)
    ind = index(trim(inputfile),'.',BACK=.TRUE.)
    filestem=inputfile(1:ind-1)

    inputformat=inputfile(ind+1:length)

    open(inputu,file=trim(inputfile),status='old',action='read')

    if(trim(inputformat)=='vasp' .or. trim(inputformat) == 'VASP') then
      call assign_str_fr_poscar(inputu,s)

    else if(trim(inputformat)=='STRUCT_OUT') then
      call assign_str_fr_struct_out(inputu,s)
    else if(trim(inputformat)=='fdf') then
      call assign_str_fr_fdf(inputu,s)
    else if(trim(inputformat)=='xsf') then
      call assign_str_fr_xsf(inputu,s)
    else
      write(*,*)
      write(*,*) 'WARNING: accepted format are vasp,arc/car,STRUCT_OUT,gjf,fdf'
      write(*,*) 'but inputformat is ',trim(inputformat)
      write(*,*) 'unrecognized input file format.'
      write(*,'(A)') 'err: check input file format.'
      stop 1
    endif
    close(inputu)

  end subroutine read_struc
end module commod
