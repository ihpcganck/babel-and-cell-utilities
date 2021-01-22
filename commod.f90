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

  subroutine show_time_now_unit(un,str)
    character(len=DSL) :: str
    integer :: un

    call date_and_time(date,time)
    WRITE(un,"(A,' date ',A4,'.',A2,'.',A2,'  ',A2,':',A2,':',A2 )") &
           trim(str),DATE(1:4),DATE(5:6),DATE(7:8),TIME(1:2),TIME(3:4),TIME(5:6)
  end subroutine show_time_now_unit

  subroutine mem_init(mem)
    type(memory_type) :: mem
    mem%n_int = 0
    mem%n_double = 0
    mem%n_dcomplex = 0
  end subroutine mem_init

  subroutine disp_tot_mem(fu,mem)
    type(memory_type) :: mem
    integer :: fu
    real(double) :: tot
    real(double),parameter :: intsize=4.0d0
    real(double),parameter :: doublesize=8.0d0
    real(double),parameter :: dcomplexsize=16.0d0
    real(double),parameter :: MegaB=1.0d6
    tot = (mem%n_int*intsize + mem%n_double*doublesize + mem%n_dcomplex*dcomplexsize)/MegaB
    write(fu,'(A,F20.10)') 'Memory (MB): ',tot
  end subroutine disp_tot_mem

  subroutine mem_inc(mem,ni,nd,nc)
    integer :: ni,nd,nc
    type(memory_type) :: mem
    mem%n_int = mem%n_int + ni
    mem%n_double = mem%n_double + nd
    mem%n_dcomplex = mem%n_dcomplex + nc
  end subroutine mem_inc

  subroutine set_up_group_op(nop,Op,gr,nfold,vhd)
    integer :: nop
    character(len=*) :: gr
    type(rmatrix33) :: Op(nop)
    integer,optional :: nfold
    character,optional :: vhd

    if(trim(gr) == 'C') then
      if(.not. present(nfold)) then
        write(*,'(A)') 'err: C is called without nfold'
        stop 1
      endif
      if(nfold <= 0) then
        write(*,'(A)') 'err: nfold is not 1 or greater for C'
        stop 1
      endif
      if(.not. present(vhd)) then
        call set_up_Cn_mat(nop,Op,nfold)
      else
        if(vhd == 'h') then
          call set_up_Cnh_mat(nop,Op,nfold)
        else if(vhd == 'v') then
          call set_up_Cnv_mat(nop,Op,nfold)
        else
          write(*,'(A)') 'err: either C_nh or C_nv.'
          stop 1
        endif
      endif
    else if(trim(gr) == 'D') then
      if(.not. present(nfold)) then
        write(*,'(A)') 'err: D is called with nfold'
        stop 1
      endif
      if(nfold <= 0) then
        write(*,'(A)') 'err: nfold is not 1 or greater for D'
        stop 1
      endif
      if(.not. present(vhd)) then
        call set_up_Dn_mat(nop,Op,nfold)
      else
        if(vhd=='h') then
          call set_up_Dh_mat(nop,Op,nfold)
        else if(vhd == 'd') then
          call set_up_Dnd_mat(nop,Op,nfold)
        else
           write(*,'(A)') 'err: There is only D_nd or D_nh'
           stop 1
        endif
      endif
    else if(trim(gr)=='Td') then
      if(nop /= 24) then
        write(*,*) 'Td: input nop is ',nop
        write(*,'(A)') 'err: nop for Td should be 24.'
        stop 1
      endif
      call set_up_Td_mat(nop,Op)
    else if(trim(gr) == 'Oh') then
      if(nop /= 48) then
        write(*,*) 'Oh: input nop is ',nop
        write(*,'(A)') 'err: nop for Oh should be 48.'
        stop 1
      endif
      call set_up_Oh_mat(nop,Op)
    else if(trim(gr) == 'Ih') then
      if(nop /= 120) then
        write(*,*) 'Ih: input nop is ',nop
        write(*,'(A)') 'err: nop for Ih should be 120.'
        stop 1
      endif
      call set_up_Ih_mat(nop,Op)
    else
      write(*,*) 'Check the group name.'
      write(*,'(A)') 'err: not implemented yet.'
      stop 1
    endif
  end subroutine set_up_group_op

  subroutine set_up_Cn_mat(nop,Op,nfold)
    integer :: nfold,nop,i
    type(rmatrix33) :: Op(nop)
    real(double) :: p(3,3)
    if(nop /= nfold) then
      write(*,*) 'nop,nfold = ',nop,nfold
      write(*,'(A)') 'err: nop and nfold should be the same.'
      stop 1
    endif
    if(nop <= 0) then
      write(*,*) 'nop in set_up_Cn_mat is ',nop
      write(*,'(A)') 'err: n must be great than 0 for C_n.'
      stop 1
    endif
    p = neta_2_rotm((/0d0,0d0,1d0/),360.0d0/(nop*one))
    Op(1)%m = p
    do i = 2, nop
      call matmatn(3,Op(i-1)%m,p,Op(i)%m)
    enddo
  end subroutine set_up_Cn_mat

  subroutine set_up_Cnh_mat(nop,Op,nfold)
    integer :: nop,nfold,i
    type(rmatrix33) :: Op(nop)
    real(double) :: xy_ref(3,3)

    if(nop /= 2*nfold) then
      write(*,*) 'nop,2*nfold=',nop,2*nfold
      write(*,'(A)') 'err: For C_nh, nop must be the same as 2*nfold'
      stop 1
    endif
    call set_up_Cn_mat(nfold,Op,nfold)
    xy_ref = reshape((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,-1d0/),(/3,3/))
    do i = 1, nfold
      call matmatn(3,Op(i)%m,xy_ref,Op(i+nfold)%m)
    enddo
  end subroutine set_up_Cnh_mat

  subroutine set_up_Cnv_mat(nop,Op,nfold)
    integer :: nop,nfold,i
    type(rmatrix33) :: Op(nop)
    real(double) :: xz_ref(3,3)
    real(double) :: M(3,3),Q(3,3),P(3,3)

    if(nop /= 2*nfold) then
      write(*,*) 'nop,2*nfold=',nop,2*nfold
      write(*,'(A)') 'err: For C_nv, nop must be the same as 2*nfold'
      stop 1
    endif
    call set_up_Cn_mat(nfold,Op,nfold)
    xz_ref = reshape((/1d0,0d0,0d0,0d0,-1d0,0d0,0d0,0d0,1d0/),(/3,3/))
    do i = 1, nfold
      call matmatn(3,Op(i)%m,xz_ref,Op(i+nfold)%m)
    enddo

    M(1:3,1:3) = reshape( (/one,zero,zero,-half,sqrt3/two,zero,zero,zero,one/), (/3,3/))

    do i = 1, 2*nfold
       write(*,*) 'i=',i
       q(1:3,1:3) = Op(i)%m(1:3,1:3)
       p = invBAB(m,q)
       write(*,*) p(1,1:3)
       write(*,*) p(2,1:3)
       write(*,*) p(3,1:3)
    enddo
  end subroutine set_up_Cnv_mat

  subroutine set_up_Dn_mat(nop,Op,nfold)
    integer :: i,nfold,nop
    type(rmatrix33) :: Op(nop)
    real(double) :: p(3,3)
    if(nop /= 2*nfold) then
      write(*,*) 'nop,2*nfold=',nop,2*nfold
      write(*,'(A)') 'err: nop and 2*nfold should be the same.'
      stop 1
    endif
    call set_up_Cn_mat(nop/2,Op,nfold)
    p = neta_2_rotm((/1.0d0,0.0d0,0.0d0/),180.0d0)
    do i = 1, nop/2
      call matmatn(3,Op(i)%m,p,Op(i+nop/2)%m)
    enddo
  end subroutine set_up_Dn_mat

  subroutine set_up_Dh_mat(nop,Op,nfold)
    integer :: nop,nfold,twonfold,i
    type(rmatrix33) :: Op(nop)
    real(double) :: xy_ref(3,3)

    if(nop /= nfold*4) then
      write(*,*) 'nop,nfold*4=',nop,nfold*4
      write(*,'(A)') 'err: for Dh, nop must be the same as 4*nfold'
      stop 1
    endif
    twonfold = 2*nfold
    call set_up_Dn_mat(twonfold,Op,nfold)
    xy_ref = reshape((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,-1d0/),(/3,3/))
    do i = 1, twonfold
      call matmatn(3,xy_ref,Op(i)%m,Op(i+twonfold)%m)
    enddo
  end subroutine set_up_Dh_mat

  subroutine set_up_Dnd_mat(nop,Op,nfold)
    integer :: nop,nfold,twonfold,i
    type(rmatrix33) :: Op(nop)
    real(double) :: p(3,3),q(3,3),r(3,3),s(3,3),t(3,3),theta
    if(nop /= nfold*4) then
      write(*,*) 'nop = ',nop
      write(*,*) 'nfold = ',nfold
      write(*,'(A)') 'err: Dnd should have 4*n elements.'
      stop 1
    endif
    twonfold = 2*nfold
    call set_up_Dn_mat(twonfold,Op,nfold)

    theta = 360d0/(2d0*nfold)/2d0
    p = neta_2_rotm((/0d0,0d0,1d0/),-theta)
    q = reshape((/1d0,0d0,0d0,0d0,-1d0,0d0,0d0,0d0,1d0/),(/3,3/))
    r = neta_2_rotm((/0d0,0d0,1d0/),theta)

    call matmatn(3,q,p,s)
    call matmatn(3,r,s,t)
    do i = 1, twonfold
      call matmatn(3,t,Op(i)%m,Op(i+twonfold)%m)
    enddo
  end subroutine set_up_Dnd_mat

  subroutine set_up_Ih_mat(nop,Op)
    integer,parameter :: ru=50
    integer :: nop,i,j
    type(rmatrix33) :: Op(nop)

    if(nop /= 120) then
      write(*,'(A)') 'err: n is not 120.'
      stop 1
    endif

    open(unit=ru,file='Ih_Op.dat',status='old',action='read')
    do i = 1, 120
      read(ru,*) j,Op(i)%m(1:3,1:3)
    enddo
    close(ru)
  end subroutine set_up_Ih_mat

  subroutine set_up_Td_mat(nop,Op)
    integer :: nop
    type(rmatrix33) :: Op(nop)
    real(double) :: x(3,3),y(3,3)

    op(1)%m = neta_2_rotm((/1.d0,1.0d0,1.0d0/),0.0d0)
    op(2)%m = neta_2_rotm((/1.d0,1.0d0,1.0d0/),120.d0)
    op(3)%m = neta_2_rotm((/1.d0,1.0d0,1.0d0/),-120.d0)
    op(4)%m = neta_2_rotm((/-1.d0,1.0d0,1.0d0/),120.d0)
    op(5)%m = neta_2_rotm((/-1.d0,1.0d0,1.0d0/),-120.d0)
    op(6)%m = neta_2_rotm((/1.d0,-1.0d0,1.0d0/),120.d0)
    op(7)%m = neta_2_rotm((/1.d0,-1.0d0,1.0d0/),-120.d0)
    op(8)%m = neta_2_rotm((/1.d0,1.0d0,-1.0d0/),120.d0)
    op(9)%m = neta_2_rotm((/1.d0,1.0d0,-1.0d0/),-120.d0)

    op(10)%m = neta_2_rotm((/1.d0,0.0d0, 0.0d0/),180.d0)
    op(11)%m = neta_2_rotm((/0.d0,1.0d0, 0.0d0/),180.d0)
    op(12)%m = neta_2_rotm((/0.d0,0.0d0, 1.0d0/),180.d0)

    x = neta_2_rotm((/1.d0,0.0d0, 0.0d0/), 90.d0)
    y = reshape((/-1,0,0,0,1,0,0,0,1/),(/3,3/))
    call matmatn(3,x,y,op(13)%m)
    x = neta_2_rotm((/1.d0,0.0d0, 0.0d0/), -90.d0)
    y = reshape((/-1,0,0,0,1,0,0,0,1/),(/3,3/))
    call matmatn(3,x,y,op(14)%m)
    x = neta_2_rotm((/0.d0,1.0d0, 0.0d0/),  90.d0)
    y = reshape((/1,0,0,0,-1,0,0,0,1/),(/3,3/))
    call matmatn(3,x,y,op(15)%m)
    x = neta_2_rotm((/0.d0,1.0d0, 0.0d0/),  -90.d0)
    y = reshape((/1,0,0,0,-1,0,0,0,1/),(/3,3/))
    call matmatn(3,x,y,op(16)%m)
    x = neta_2_rotm((/0.d0,0.0d0, 1.0d0/),   90.d0)
    y = reshape((/1,0,0,0, 1,0,0,0,-1/),(/3,3/))
    call matmatn(3,x,y,op(17)%m)
    x = neta_2_rotm((/0.d0,0.0d0, 1.0d0/),  -90.d0)
    y = reshape((/1,0,0,0, 1,0,0,0,-1/),(/3,3/))
    call matmatn(3,x,y,op(18)%m)

    op(19)%m = reshape((/0,1,0,1,0,0,0,0,1/),(/3,3/))
    op(20)%m = reshape((/0,-1,0,-1,0,0,0,0,1/),(/3,3/))
    op(21)%m = reshape((/1, 0,0, 0,0,1,0,1,0/),(/3,3/))
    op(22)%m = reshape((/1, 0,0, 0,0,-1,0,-1,0/),(/3,3/))
    op(23)%m = reshape((/0, 0,1, 0,1, 0,1, 0,0/),(/3,3/))
    op(24)%m = reshape((/0, 0,-1, 0,1, 0,-1, 0,0/),(/3,3/))

  end subroutine set_up_Td_mat

  subroutine set_up_Oh_mat(nop,Op)
    integer :: nop,i
    type(rmatrix33) :: Op(nop)
    real(double) :: Inversion(3,3)

    if(nop /= 48) then
      write(*,*) 'nop = ',nop
      write(*,'(A)') 'err: nop should be 48 for Oh.'
      stop 1
    endif
    call set_up_Td_mat(24,Op)
    Inversion = reshape((/-1,0,0,0,-1,0,0,0,-1/),(/3,3/))
    do i = 1, 24
      call matmatn(3,Inversion,op(i)%m,op(i+24)%m)
    enddo
  end subroutine set_up_Oh_mat

  subroutine inversemapping(n,fdindex,bdindex)
    integer :: n,i,ind
    integer :: fdindex(n)
    integer :: bdindex(n)
    do i = 1, n
      ind = fdindex(i)
      bdindex(ind) = i
    enddo
  end subroutine inversemapping

  function num2str(v,nspaces) result (str)
    integer :: tmpv,i,v,nspaces,digit,ind
    character(len=nspaces) :: str
    tmpv = v
    if(v < 0) then
      write(*,'(A)') 'err: v is negative.'
      stop 1
    endif
    do i = 1, nspaces
      digit = modulo(tmpv,10)
      tmpv = tmpv/10
      ind = nspaces-i+1
      str(ind:ind) = DIGITARR(digit)
    enddo
  end function num2str

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

  subroutine lowercase(string)
    integer :: i,j
    character(len=*)  :: string
    do i = 1,len(string)
      j = index(upper,string(i:i))
      if(j /= 0) string(i:i) = lower(j:j)
    enddo
  end subroutine lowercase

  function num_of_lines(funit,filename)
    integer :: funit
    character(len=*) :: filename
    integer :: num_of_lines
    integer :: nlines,io

    nlines = 0
    open(funit,file=trim(filename),status='old',action='read')
    do

      read(funit,*,iostat=io)
      if(io /= 0) exit
      nlines = nlines + 1
    enddo
    close(funit)
    num_of_lines = nlines
  end function num_of_lines

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

  function newunit()
    integer :: newunit
    integer, parameter :: lower_unitn=10, upper_unitn=1000
    logical :: is_opened
    integer :: i

    newunit = -1
    do i = lower_unitn, upper_unitn
      inquire(unit=i,opened=is_opened)
      if (.not. is_opened) then
        newunit = i
        exit
      endif
    enddo
    if(newunit == -1) then
      write(*,*) 'lower_unitn,upper_unitn=',lower_unitn,upper_unitn
      write(*,'(A)') 'err: cannot find a unused unit number.'
      stop 1
    endif
  end function newunit

  subroutine seq_read_gendata(dstype,fileu,nametag,v)
    integer :: dstype,fileu
    character(len=*) :: nametag
    type(dtype) :: v
    character(len=DSL) :: tmpstr

    write(*,*)
    write(*,'(A)') 'Looking for keyword = '//trim(nametag)

    v%li = .false.
    v%lr = .false.
    v%ls = .false.
    if(dstype == itype) then
      v%li = .true.
    else if(dstype == rtype) then
      v%lr = .true.
    else if(dstype == stype) then
      v%ls = .true.
    endif

    if(v%li) then
      read(fileu,*) tmpstr,v%i
    else if(v%lr) then
      read(fileu,*) tmpstr,v%r
    else if(v%ls) then
      read(fileu,*) tmpstr,v%s
    endif
    if(trim(tmpstr) /= trim(nametag)) then
      write(*,*) 'tmpstr is '//trim(tmpstr)
      write(*,*) 'keyword '//trim(nametag)//' is expected.'
      write(*,'(A)') 'err: seq_read_gendata. check keyword.'
      stop 1
    endif
    if(v%li) then
      write(*,'(A,I10)') 'Integer type found. The assigned variable is ',v%i
    else if(v%lr) then
      write(*,'(A,E20.13)') 'Real type found. The assigned variable is ',v%r
    else if(v%ls) then
      write(*,'(A)') 'String type found. The assigned variable is '//trim(v%s)
    endif
  end subroutine seq_read_gendata

  subroutine read_gendata(dstype,fileu,filename,nametag,v,dopt,vdefault,vec)
    integer,optional :: vec
    integer :: dopt
    integer :: dstype,fileu
    type(dtype) :: v,vdefault
    character(len=*) :: nametag,filename
    logical :: done
    integer :: ios,s
    character(len=DSL) :: str1,str2
    logical :: found
    integer :: i,j
    character(len=DSL) :: tmpstr

    if(present(vec)) then
      if(vec > DTYPEMAXLEN) then
        write(*,*) 'vec= ',vec
        write(*,*) 'DTYPEMAXLEN= ',DTYPEMAXLEN
        write(*,'(A)') 'err: vec > DTYPEMAXLEN.'
        stop 1
      endif
    endif

    found = .false.
    v%li = .false.
    v%lr = .false.
    v%ls = .false.
    if(dstype == itype) then
      v%li = .true.
    else if(dstype == rtype) then
      v%lr = .true.
    else if(dstype == stype) then
      v%ls = .true.
    endif

    write(*,'(A)') 'Search keyword=|'//trim(nametag)//'| in '//trim(filename)
    open(fileu,file=trim(filename),status='old',action='read')

    done = .false.
    do
      if(done) exit
      read(fileu,DSF,iostat=ios) str1
      if(ios /= 0) then
        done = .true.
      else
        call onespace(DSL,str1,str2)

        s = scan(str2,' ')
        if(s > 0) then
          if(str2(1:s-1) == '#keyword') then
            str2 = str2(s+1:DSL)
            s = scan(str2,' ')
            if(str2(1:s-1) == trim(nametag)) then
              done = .true.
              found = .true.
              str2 = str2(s+1:DSL)
              if(v%li) then

                write(*,'(A)') 'str2 is |'//trim(str2)//'|'
                if(present(vec)) then
                  read(str2,*) v%ivec(1:vec)
                  tmpstr = N2str(v%ivec(1))
                  do j = 2, vec
                    tmpstr = trim(tmpstr)//','//trim(N2str(v%ivec(j)))
                  enddo
                  write(*,'(A)') 'Integer variable '//trim(nametag)//' is |'//trim(tmpstr)//'|'
                else
                  read(str2,*) v%i
                  write(*,'(A)') 'Integer variable '//trim(nametag)//' is '//trim(N2str(v%i))
                endif
              else if(v%lr) then

                write(*,'(A)') 'str2 is |'//trim(str2)//'|'
                if(present(vec)) then
                  read(str2,*) v%rvec(1:vec)
                  write(*,*) 'Real variable '//trim(nametag)//' is ',v%rvec(1:vec)
                else
                  read(str2,*) v%r
                  write(*,*) 'Real variable '//trim(nametag)//' is ',v%r
                endif
              else if(v%ls) then

                write(*,'(A)') 'str2 is |'//trim(str2)//'|'
                if(present(vec)) then
                  read(str2,*) v%svec(1:vec)
                  write(*,*) 'String variable '//trim(nametag)//' is:'
                  do i = 1, vec
                    write(*,'(A)') trim(v%svec(i))
                  enddo
                else
                  read(str2,*) v%s
                  write(*,'(A)') 'String variable '//trim(nametag)//' is '//trim(v%s)
                endif
              endif
            endif
          else

          endif
        else
          done = .true.
        endif
      endif
    enddo
    close(fileu)

    if(.not. found .and. dopt == NONDEFAULTOPT) then
      write(*,*) 'cannot find the keyword='//trim(nametag)
      write(*,'(A)') 'err: cannot find the the keyword.'
      stop 1
    else if(.not. found .and. dopt == DEFAULTOPT) then
      write(*,'(/A)') 'Could not find '//trim(nametag)
      if(v%li) then
        if(present(vec)) then
          v%ivec(1:vec) = vdefault%ivec(1:vec)
          write(*,*) 'Warning: assume default value of ',v%ivec(1:vec)
        else
          v%i = vdefault%i
          write(*,*) 'Warning: assume default value of ',v%i
        endif
      else if(v%lr) then
        if(present(vec)) then
          v%rvec(1:vec) = vdefault%rvec(1:vec)
          write(*,*) 'Warning: assume default value of ',v%rvec(1:vec)
        else
          v%r = vdefault%r
          write(*,*) 'Warning: assume default value of ', v%r
        endif
      else if(v%ls) then
        if(present(vec)) then
          v%svec(1:vec) = vdefault%svec(1:vec)
          write(*,*) 'Warning: assume default value of ',v%svec(1:vec)
        else
          v%s = vdefault%s
          write(*,'(A)') 'Warning: assume default value of '//trim(v%s)
        endif
      endif
    endif
  end subroutine read_gendata

  function sppolar2cart(p) result(r)
    real(double) :: rmag,theta,phi,p(3),r(3)
    rmag = p(1)
    theta = p(2)*deg2rad
    phi = p(3)*deg2rad
    r(1:3) = (/rmag*sin(theta)*cos(phi),rmag*sin(theta)*sin(phi),rmag*cos(theta)/)
  end function sppolar2cart

  function rthetaphi(p)
    real(double) :: p(3),rthetaphi(3),testv(3),diffv(3)
    real(double) :: rmag,phi,theta,rho,q(3),pmax,diff

    real(double),parameter :: mag_tol=1.0d-12

    real(double),parameter :: eps=1.0d-11

    q(1:3) = p(1:3)
    rmag = vecmag3(p)
    if(rmag < mag_tol) then
      write(*,*) 'p = ',p
      write(*,*) 'rmag = ',rmag
      write(*,*) 'mag_tol = ',mag_tol
      write(*,'(A)') 'err: The magnitude of p is too small (less than mag_tol)'
      stop 1
    else

      pmax = abs(p(1))
      if(abs(p(2)) > pmax) then
        pmax = abs(p(2))
      endif
      if(abs(p(3)) > pmax) then
        pmax = abs(p(3))
      endif

      q(1:3) = q(1:3)/pmax

    endif

    phi = atan2(q(2),q(1))

    if(sqrt(q(1)*q(1) + q(2)*q(2)) < 1.0d-12) then

      phi = zero

    endif

    if(phi < 0) then
      phi = 2*pi + phi
    endif

    rho = sqrt(q(1)*q(1)+q(2)*q(2))
    theta = atan2(rho,q(3))

    rthetaphi(1:3) = (/rmag,theta,phi/)

    testv(1:3) = (/rmag*sin(theta)*cos(phi),rmag*sin(theta)*sin(phi),rmag*cos(theta)/)
    diffv(1:3) = testv(1:3) - p(1:3)
    diff = vecmag3(diffv)
    if(diff > eps) then
      write(*,*) 'eps,diff= ',eps,diff
    endif

  end function rthetaphi

  function neta_2_rotm(n,eta) result (m)
    real(double) :: m(3,3),n(3),eta,a,b,c,ct,st,nmag

    ct = cos(eta*deg2rad)
    st = sin(eta*deg2rad)

    nmag = sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
    if(abs(nmag) < 1.0d-5) then
      write(*,'(A)') 'err: in rotmat: mag u is too small.'
      stop 1
    endif
    a = n(1)/nmag; b = n(2)/nmag; c = n(3)/nmag
    m(1,1) = a*a*(1-ct)+ct; m(1,2) = a*b*(1-ct)-c*st; m(1,3) = a*c*(1-ct)+b*st
    m(2,1) = a*b*(1-ct)+c*st; m(2,2) = b*b*(1-ct)+ct; m(2,3) = b*c*(1-ct)-a*st
    m(3,1) = a*c*(1-ct)-b*st; m(3,2) = b*c*(1-ct)+a*st; m(3,3) = c*c*(1-ct)+ct
  end function neta_2_rotm

  function euler_2_rotm(alpha,beta,gamm) result(m)
    real(double) :: alpha,beta,gamm
    real(double) :: a,b,g
    real(double) :: m(3,3)
    real(double) :: sa,ca,sb,cb,sg,cg

    a = alpha*deg2rad
    b = beta*deg2rad
    g = gamm*deg2rad
    sa = sin(a)
    ca = cos(a)
    sb = sin(b)
    cb = cos(b)
    sg = sin(g)
    cg = cos(g)
    m(1,1:3) = (/ca*cb*cg-sa*sg,-ca*cb*sg-sa*cg,ca*sb/)
    m(2,1:3) = (/sa*cb*cg+ca*sg,-sa*cb*sg+ca*cg,sa*sb/)
    m(3,1:3) = (/-sb*cg,sb*sg,cb/)
  end function euler_2_rotm

  function euler_2_neta(EulerAngle) result (p)
    real(double) :: EulerAngle(3),p(4)
    real(double) :: heta,alpha,beta,gamm,cosheta,sinheta,eta,a,b,c

    alpha = EulerAngle(1)*deg2rad
    beta = EulerAngle(2)*deg2rad
    gamm = EulerAngle(3)*deg2rad
    cosheta = cos(beta/two)*cos((alpha+gamm)/two)

    sinheta = sqrt(one - cosheta**two)
    heta = atan2(sinheta,cosheta)
    eta = two*heta

    if(eta < 1.0d-12) then
      a = one; b = one; c = one
    else

      a = - sin(beta/two)*sin((alpha-gamm)/two)/sin(eta/two)
      b = sin(beta/two)*cos((alpha-gamm)/two)/sin(eta/two)
      c = cos(beta/two)*sin((alpha+gamm)/two)/sin(eta/two)
    endif

    p(1) = a
    p(2) = b
    p(3) = c
    p(4) = eta*rad2deg
  end function euler_2_neta

  function rotm_2_euler(m) result(EulerAngle)
    real(double) :: sinbeta,m(3,3),EulerAngle(3)
    real(double) :: alpha,beta,gamm,d1pd4,d2md3,sum_alpha_gamm
    real(double),parameter :: eps=1.0d-12
    real(double) :: diff,diffmat(3,3),checkmat(3,3)
    integer :: toprint

    sinbeta = sqrt(m(1,3)**2 + m(2,3)**2)

    gamm = atan2(m(3,2),-m(3,1))*rad2deg
    alpha = atan2(m(2,3),m(1,3))*rad2deg
    beta = atan2(sinbeta,m(3,3))*rad2deg

    toprint = 0

    if(toprint == 1) then
      write(*,*) 'sinbeta,s(3,3)= ',sinbeta,m(3,3)
      write(*,*) 'alpha = ',alpha
      write(*,*) 'beta = ',beta
      write(*,*) 'gamm = ', gamm

      if(abs(sinbeta) < eps) then
        write(*,*) 'sinbeta= ',sinbeta
        write(*,*) 'WARNING: sinbeta is too small. But we are okay if we do not handle near 2-fold rotations'
      endif
    endif

    if(abs(beta) > 1.0d-10) then

    else
      if(toprint == 1) then
        write(*,*) 'WARNING: beta is too close to zero, we symmetrize alpha and gamma'
        write(*,*) 'before special case handling...'
        write(*,*) 'alpha = ',alpha
        write(*,*) 'beta = ',beta
        write(*,*) 'gamm = ', gamm
      endif

      d1pd4 = m(1,1)+m(2,2)
      d2md3 = m(2,1)-m(1,2)
      sum_alpha_gamm = atan2(d2md3,d1pd4)
      alpha = (sum_alpha_gamm/2d0)*rad2deg
      gamm = (sum_alpha_gamm/2d0)*rad2deg

      toprint = 0
      if(toprint == 1) then
        write(*,*) 'after special case handling...'
        write(*,*) 'rotm_2_euler: alpha = ',alpha
        write(*,*) 'rotm_2_euler: beta = ',beta
        write(*,*) 'rotm_2_euler: gamm = ', gamm
      endif
    endif

    EulerAngle(1:3) = (/alpha,beta,gamm/)

    checkmat(1:3,1:3) = euler_2_rotm(EulerAngle(1),EulerAngle(2),EulerAngle(3))

    diffmat(1:3,1:3) = checkmat(1:3,1:3) - m(1:3,1:3)
    diff = matmagn(3,diffmat(1,1))

    if(diff > eps) then
      write(*,*) 'rotm_2_euler: eps,diff= ',eps,diff
    endif
  end function rotm_2_euler

  function neta_2_euler(p) result (EulerAngle)
    real(double) :: n(3),eta,p(4),EulerAngle(3),m(3,3)
    n(1:3) = p(1:3)
    eta = p(4)
    m = neta_2_rotm(n,eta)
    EulerAngle = rotm_2_euler(m)
  end function neta_2_euler

  function rotm_2_neta(m) result(p)
    real(double) :: m(3,3),p(4),euler(3), checkm(3,3),diffm(3,3)
    real(double) :: eps
    integer :: info
    integer,parameter :: ldvl=3
    integer,parameter :: ldvr=3
    complex(double) :: VL(ldvl,3),VR(ldvr,3)
    integer,parameter :: LCWORK=6
    integer,parameter :: LRWORK=6
    complex(double) :: CWORK(LCWORK),cmat(3,3)
    real(double) :: RWORK(LRWORK)
    complex(double) :: eigval(3),sumeigval
    integer :: i
    real(double) :: rsum,costheta,sintheta,theta,xsq,tracem
    logical :: found,near2fold
    real(double) :: diff
    integer :: toprint

    eps = 1.0d-7

    call check_orthogonal(m(1,1),eps)

    tracem = m(1,1) + m(2,2) + m(3,3)

    cmat(1:3,1:3) = m(1:3,1:3)

    call zgeev('V','V',3,cmat(1,1),3,eigval(1),VL,LDVL,VR,LDVR,CWORK(1),LCWORK,RWORK(1),info)
    if(info /= 0) then
      write(*,'(A)') 'err: error in zgeev.'
      stop 1
    endif

    sumeigval = eigval(1)+eigval(2)+eigval(3)
    rsum = real(sumeigval)

    if(abs(rsum - tracem) > eps) then
      write(*,*) 'rsum,tracem= ',rsum,tracem
      write(*,'(A)') 'err: rsum and tracem must be the same.'
      stop 1
    endif

    if(abs ( abs(tracem) - 3.0d0)  < 1.0d-12) then
      write(*,*) 'Hit near identity: return a reasonabe n and eta'
      p(1:3) = one
      p(4) = zero
      return
    endif

    near2fold = .false.

    if(abs(rsum - (-one)) < eps) then

      found = .false.
      do i = 1, 3

        if(found) cycle

        if( abs (real(eigval(i) - one))  < eps) then
          found = .true.
          near2fold = .true.

          costheta = (rsum- one)/two
          xsq = one- costheta*costheta
          if(xsq < zero) then

            xsq = -xsq
          endif
          sintheta = sqrt(xsq)
          theta = atan2(sintheta,costheta)
          p(4) = theta*rad2deg
          p(1) = real(vr(1,i))
          p(2) = real(vr(2,i))
          p(3) = real(vr(3,i))
        endif
      enddo
    endif

    if(.not. near2fold) then

      euler = rotm_2_euler(m)

      p = euler_2_neta(euler)
    endif

    checkm = neta_2_rotm( (/p(1),p(2),p(3)/), p(4) )
    diffm = checkm-m
    diff = matmagn(3,diffm(1,1))

    toprint = 0
    if(toprint == 1) then
      write(*,*) 'm,checkm,diffm='
      write(*,'(3F8.4,A,3F8.4,A,3F8.4)') m(1,1:3), '  |  ', checkm(1,1:3),'  |  ',diffm(1,1:3)
      write(*,'(3F8.4,A,3F8.4,A,3F8.4)') m(2,1:3), '  |  ', checkm(2,1:3),'  |  ',diffm(2,1:3)
      write(*,'(3F8.4,A,3F8.4,A,3F8.4)') m(3,1:3), '  |  ', checkm(3,1:3),'  |  ',diffm(3,1:3)
      write(*,*) 'rotm_2_neta: eps,diff= ',eps,diff
    endif
    if(diff > eps) then
      write(*,'(A)') 'err: in rotm_2_neta: m and p not consistent.'
      stop 1
    endif
  end function rotm_2_neta

  subroutine GetLenAng(a,L)
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
  end subroutine GetLenAng

  subroutine getreallatt(L,a)
    real(double) :: a(3,3),L(6),ang(3)
    real(double) :: cz(3,3)
    integer :: i
    integer,parameter :: PRINT_ALTERNATIVE=100
    integer,parameter :: noPRINT_ALTERNATIVE=101
    integer :: print_option

    do i = 1, 3
      if(L(i) <= 0) then
        write(*,*) 'getreallatt: i=',i,' ,L(i) = ',L(i)
        write(*,'(A)') 'err: Check L(i) value.'
        stop 1
      endif
    enddo

    do i = 1, 3
      ang(i) = L(i+3)

      if(ang(i) <= 0.0d0 .or. ang(i) >= 180.0d0) then
        write(*,*) 'GenReaLatt, i,ang(i)=',i,ang(i)
        write(*,'(A)') 'err: angle problem.'
        stop 1
      endif
      ang(i) = ang(i)*deg2rad
    enddo
    a(:,:) = 0.0D0
    a(1,1) = L(1)
    a(1,2) = L(2)*cos(ang(3))
    a(2,2) = L(2)*sin(ang(3))
    a(1,3) = L(3)*cos(ang(2))
    a(2,3) = (L(2)*L(3)*cos(ang(1)) - a(1,2)*a(1,3))/a(2,2)
    a(3,3) = sqrt(L(3)*L(3)-a(1,3)*a(1,3)-a(2,3)*a(2,3))

    print_option=noPRINT_ALTERNATIVE

    if(print_option == PRINT_ALTERNATIVE) then

      cz(:,:) = zero
      cz(1:3,1) = (/ L(1)*sin(ang(2)), zero, L(1)*cos(ang(2)) /)
      cz(1:3,2) = L(2)* (/ ( cos(ang(3)) - cos(ang(1))*cos(ang(2)))/sin(ang(2)),  &
                  sqrt( one  - cos(ang(1))**2 - cos(ang(2))**2 - cos(ang(3))**2 + &
                  two* cos(ang(1))*cos(ang(2)) * cos(ang(3)))/sin(ang(2)),  cos(ang(1)) /)
      cz(1:3,3) = (/zero, zero, L(3) /)
      write(*,'(A,3F20.10)') 'alternative: a1 is ',cz(1:3,1)
      write(*,'(A,3F20.10)') 'alternative: a2 is ',cz(1:3,2)
      write(*,'(A,3F20.10)') 'alternative: a3 is ',cz(1:3,3)
    endif

  end subroutine getreallatt

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

  function angle_between_vecs3(A,B) result (C)
    real(double) :: A(3),B(3),C,maga,magb,AdotB,xsq,costheta,sintheta
    real(double),parameter :: eps=1.0d-8

    maga = vecmag3(A(1))
    magb = vecmag3(B(1))

    if(maga < eps .or. magb < eps) then
      write(*,'(A)') 'err: maga or magb too small.'
      stop 1
    endif
    adotb = dotprod3(A(1),B(1))
    costheta = adotb/(maga*magb)
    xsq = one-costheta*costheta
    if(xsq < zero) then
      xsq = -xsq
    endif
    sintheta = sqrt(xsq)
    c = atan2(sintheta,costheta)
    c = c*rad2deg
  end function angle_between_vecs3

  function CrossProd(A,B)  result(C)
    real(double) :: A(1:3),B(1:3),C(1:3)
    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1) - A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)
  end function CrossProd

  function VecProd(A,B) result(C)
    real(double) :: A(1:3),B(1:3),C(1:3)
    C = CrossProd(A,B)
  end function VecProd

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

  function TripProd(A)
    real(double) :: A(3,3),TripProd
    TripProd = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - &
               A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + &
               A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  end function TripProd

  function det3(a)
    real(double) :: det3,a(3,3)
    det3 = tripprod(a)
  end function det3

  function tracen(n,a)
    integer :: n,i
    real(double) :: tracen,a(n,n)
    tracen = zero
    do i = 1, n
      tracen = tracen + a(i,i)
    enddo
  end function tracen

  function trace3(a)
    real(double) :: trace3,a(3,3)
    trace3 = tracen(3,a(1,1))
  end function trace3

  function ITripProd(A)
    integer :: A(3,3),ITripProd
    ITripProd = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - &
               A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + &
               A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  end function ITripProd

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
    call getLenAng(s2%a,s2%la)
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

  subroutine put_in_bbox(pd,boxpd,xyzbound)
    type(supercell) :: pd,boxpd
    real(double) :: xyzbound(3,2)
    real(double) :: p(3),xgap,ygap,zgap
    integer :: i,j

    xyzbound(1:3,1) = 1d100; xyzbound(1:3,2) = -1d100
    do i = 1, pd%n
      do j = 1, 3
        if(pd%at(i)%ac(j) < xyzbound(j,1)) then
          xyzbound(j,1) = pd%at(i)%ac(j)
        endif
        if(pd%at(i)%ac(j) > xyzbound(j,2)) then
          xyzbound(j,2) = pd%at(i)%ac(j)
        endif
      enddo
    enddo

    xgap = 5d0; ygap = 5d0; zgap = 5d0

    boxpd%la = (/xyzbound(1,2)-xyzbound(1,1)+xgap,xyzbound(2,2)-xyzbound(2,1)+ygap,xyzbound(3,2)-xyzbound(3,1)+zgap,90d0,90d0,90d0/)
    call getreallatt(boxpd%la,boxpd%a)
    call real2recip(boxpd%a,boxpd%b)
    boxpd%n = pd%n
    do i = 1, pd%n
      p = pd%at(i)%ac - (/xyzbound(1,1),xyzbound(2,1),xyzbound(3,1)/)
      p = p + 0.5*(/xgap,ygap,zgap/)
      call abs2frac(boxpd%b,p,boxpd%at(i)%f)
      boxpd%at(i)%z = pd%at(i)%z
    enddo
  end subroutine put_in_bbox

  subroutine epsilondiff(a,b,eps)
    real(double) :: a,b,eps
    real(double) :: cp,cm
    if(abs(a) >= abs(b)) then
      cp = a
      cm = b
    else
      cp = b
      cm = a
    endif
    if(cp >= 0) then
      if(cm >= cp*(1d0-eps) .and. cm <= cp*(1d0+eps)) then

      else
        write(*,*) 'cp,cm,eps=',cp,cm,eps
        write(*,*) 'cp*(1d0-eps) = ',cp*(1d0-eps)
        write(*,*) 'cp*(1d0+eps) = ',cp*(1d0+eps)
        write(*,'(A)') 'err: cp is in the non-negative loop, a and b are different.'
        stop 1
      endif
    else
      if(cm >= cp*(1d0+eps) .and. cm <= cp*(1d0-eps)) then

      else
        write(*,*) 'cp,cm,eps=',cp,cm,eps
        write(*,*) 'cp*(1d0+eps) = ',cp*(1d0+eps)
        write(*,*) 'cp*(1d0-eps) = ',cp*(1d0-eps)
        write(*,'(A)') 'err: cp is in the negative loop, a and b are different.'
        stop 1
      endif
    endif
  end subroutine epsilondiff

  function matmagn(n,m)
    integer :: n
    real(double) :: matmagn,m(n,n)
    integer :: i,j
    matmagn = zero
    do j = 1, n
      do i = 1, n
        matmagn = matmagn + m(i,j)**2
      enddo
    enddo
    matmagn = sqrt(matmagn)
  end function matmagn

  function zmatmagn(n,m)
    integer :: i,j,n
    complex(double) :: m(n,n)
    real(double) :: matmag,zmatmagn

    matmag = zero
    do j = 1, n
      do i = 1, n
        matmag = matmag + real(m(i,j)*conjg(m(i,j)))
      enddo
    enddo
    matmag = sqrt(matmag)
    zmatmagn = matmag
  end function zmatmagn

  function Frobenius3(A)
    real(double) :: Frobenius3,A(3,3),elesum
    integer :: i,j
    elesum = zero
    do i = 1, 3
      do j = 1, 3
        elesum = elesum + A(j,i)**2
      enddo
    enddo
    Frobenius3 = sqrt(elesum)
  end function Frobenius3

  function vecmag3(x)
    real(double) :: vecmag3,x(3)
    vecmag3 = sqrt(dotprod3(x,x))
  end function vecmag3

  function vecmag3_pytha(x)
    real(double) :: a12,vecmag3_pytha,x(3)
    a12 = pythagoras(x(1),x(2))
    vecmag3_pytha = pythagoras(a12,x(3))
  end function vecmag3_pytha

  function pythagoras(a,b)
    real(double) :: a,b,pythagoras
    real(double) :: p,q,r,s,t
    integer :: iter
    p = max(abs(a),abs(b))
    q = min(abs(a),abs(b))
    if(abs(q) == 0d0) then
      pythagoras = p
    else
      iter = 0
      do
        iter = iter + 1
        r = (q/p)**2
        t = 4d0 + R
        if(t  == 4d0) exit
        s = r/t
        p = p + 2d0*p*s
        q = q*s
      enddo
      pythagoras = p
      write(*,*) 'iter = ',iter
    endif
  end function pythagoras

  function vecmagn(n,x)
    integer :: n,i
    real(double) :: mag,vecmagn,x(n)
    mag = 0d0
    do i = 1, n
      mag = mag + x(i)**2
    enddo
    vecmagn = sqrt(mag)
  end function vecmagn

  subroutine LDLT(n,A,L,D,v)
    integer :: n,i,j,k
    real(double) :: A(n,n),L(n,n),D(n,n),v(n)
    real(double) :: sumj,sumk
    do i = 1, n
      do j = 1, n
        L(j,i) = zero
        D(j,i) = zero
      enddo
    enddo
    do i = 1, n
      L(i,i) = 1
    enddo
    do i = 1, n
      do j = 1, i-1
        v(j) = L(i,j)*D(j,j)
      enddo
      sumj = zero
      do j = 1, i-1
        sumj = sumj + L(i,j)*v(j)
      enddo
      d(i,i) = A(i,i) - sumj
      do j = i+1, n
        sumk = zero
        do k = 1, i-1
          sumk = sumk + L(j,k)*v(k)
        enddo
        L(j,i) = (A(j,i) - sumk)/d(i,i)
      enddo
    enddo
  end subroutine LDLT

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

  function inv2x2f(a) result(ainv)
    real(double) :: a(2,2),ainv(2,2)
    call inv2x2(a(1,1),ainv(1,1))
  end function inv2x2f

  subroutine inv2x2(p,q)
    real(double) :: p(2,2),q(2,2),a,b,c,d,det
    a = p(1,1)
    b = p(1,2)
    c = p(2,1)
    d = p(2,2)
    q(1,1) = d
    q(2,2) = a
    q(1,2) = -b
    q(2,1) = -c
    det = a*d - b*c
    if(abs(det) < 1d-10) then
      write(*,'(A)') 'err: 2x2 matrix inversion, det is near zero.'
      stop 1
    endif
    q = q/det
  end subroutine inv2x2

  subroutine complexHCn(n,A,B)
    integer :: i,j,n
    complex(double) :: A(n,n),B(n,n)
    do i = 1, n
      do j = 1, n
        B(i,j) = conjg(A(j,i))
      enddo
    enddo
  end subroutine complexHCn

  subroutine transposen(n,A,B)
    integer :: i,j,n
    real(double) :: A(n,n),B(n,n)
    do i = 1, n
      do j = 1, n
        B(i,j) = A(j,i)
      enddo
    enddo
  end subroutine transposen

  function transp3(A) result (AT)
    real(double) :: A(3,3), AT(3,3)
    AT = transpn(3,A(1,1))
  end function transp3

  function transpn(n,A) result (AT)
    integer :: n
    real(double) :: A(n,n), AT(n,n)
    call transposen(n,A(1,1),AT(1,1))
  end function transpn

  subroutine matpmatn(n,a,b,c)
    integer :: i,j,n
    real(double) :: a(n,n),b(n,n),c(n,n)
    do i = 1, n
      do j = 1, n
          c(j,i) = a(j,i)+b(j,i)
      enddo
    enddo
  end subroutine matpmatn

  subroutine matmat3(a,b,c)
    real(double):: a(3,3),b(3,3),c(3,3)
    call matmatn(3,a(1,1),b(1,1),c(1,1))
  end subroutine matmat3

  subroutine matmatn(n,a,b,c)
    integer :: i,j,k,n
    real(double) :: a(n,n),b(n,n),c(n,n),sumr
    do i = 1, n
      do j = 1, n
        sumr = zero
        do k = 1,n
          sumr = sumr + a(i,k)*b(k,j)
        enddo
        c(i,j) = sumr
      enddo
    enddo
  end subroutine matmatn

  subroutine matmatmnp(m,n,p,a,b,c)
    integer :: i,j,k,n,m,p
    real(double) :: a(m,n),b(n,p),c(m,p)

    do i = 1, m
      do j = 1, p
        c(i,j) = 0D0
        do k = 1,n
          c(i,j) = c(i,j) + a(i,k)*b(k,j)
        enddo
      enddo
    enddo
  end subroutine matmatmnp

  subroutine matpmatn_update(n,a,b)
    integer :: i,j,n
    real(double) :: a(n,n),b(n,n)
    do i = 1, n
      do j = 1, n
          a(j,i) = a(j,i)+b(j,i)
      enddo
    enddo
  end subroutine matpmatn_update

  subroutine zmatmatn(n,a,b,c)
    integer :: i,j,k,n
    complex(double) :: a(n,n),b(n,n),c(n,n)
    do i = 1, n
      do j = 1, n
        c(i,j) = 0D0
        do k = 1,n
          c(i,j) = c(i,j) + a(i,k)*b(k,j)
        enddo
      enddo
    enddo
  end subroutine zmatmatn

  subroutine zmatvecn(n,a,b,c)
    integer :: j,k,n
    complex(double) :: a(n,n),b(n),c(n)
    do j = 1, n
      c(j) = czero
      do k = 1,n
        c(j) = c(j) + a(j,k)*b(k)
      enddo
    enddo
  end subroutine zmatvecn

  subroutine matvecn(n,a,b,c)
    integer :: j,k,n
    real(double) :: a(n,n),b(n),c(n)
    do j = 1, n
      c(j) = 0D0
      do k = 1,n
        c(j) = c(j) + a(j,k)*b(k)
      enddo
    enddo
  end subroutine matvecn

  subroutine matvec3(a,b,c)
    real(double) :: a(3,3),b(3),c(3)
    call matvecn(3,a(1,1),b(1),c(1))
  end subroutine matvec3

  subroutine matvecmn(m,n,a,b,c)
    integer :: j,k,m,n
    real(double) :: a(m,n),b(n),c(m)
    do j = 1, m
      c(j) = 0D0
      do k = 1,n
        c(j) = c(j) + a(j,k)*b(k)
      enddo
    enddo
  end subroutine matvecmn

  function rotmat_x(t) result(m)
    real(double) :: m(3,3),t,u(3)
    u(1:3) = (/1,0,0/)
    m = neta_2_rotm(u,t)
  end function rotmat_x

  function rotmat_y(t) result(m)
    real(double) :: m(3,3),t,u(3)
    u(1:3) = (/0,1,0/)
    m = neta_2_rotm(u,t)
  end function rotmat_y

  function rotmat_z(t) result(m)
    real(double) :: m(3,3),t,u(3)
    u(1:3) = (/0,0,1/)
    m = neta_2_rotm(u,t)
  end function rotmat_z

  subroutine check_orthogonal(R,eps)
    real(double) :: diffabs,R(3,3),eps,tmp1(3,3),tmp2(3,3),tmp3(3,3),diffm(3,3)
    integer :: i,j

    do i = 1, 3
      do j = 1, 3
        tmp1(i,j) = R(j,i)
      enddo
    enddo
    call matmat3(R(1,1),tmp1(1,1),tmp2(1,1))
    tmp3 = zero
    tmp3(1,1) = one
    tmp3(2,2) = one
    tmp3(3,3) = one
    diffm = tmp2 - tmp3
    diffabs = matmagn(3,diffm(1,1))
    if(diffabs > eps) then
      write(*,*) ' R= '
      write(*,'(3F12.5)') R(1,1:3)
      write(*,'(3F12.5)') R(2,1:3)
      write(*,'(3F12.5)') R(3,1:3)
      write(*,*) 'diffabs = ',diffabs
      write(*,*) 'Are you sure that the matrix R is orthogonal ?'
      write(*,'(A)') 'err: R is not orthogonal within eps.'
      stop 1
    endif
  end subroutine check_orthogonal

  function invBAB(B,A) result (m)
    real(double) :: B(3,3),A(3,3),m(3,3),tmp(3,3)
    call matmatn(3,inv3x3(B),A,tmp)
    call matmatn(3,tmp,B,m)
  end function invBAB

  function BAinvB(B,A) result (m)
    real(double) :: B(3,3),A(3,3),m(3,3),tmp(3,3)
    call matmatn(3,B,A,tmp)
    call matmatn(3,tmp,inv3x3(B),m)
  end function BAinvB

  function invRAR(R,A) result (m)
    real(double) :: invr(3,3),R(3,3),A(3,3),m(3,3)

    invr = inv3x3(r(1,1))
    m = RAinvR(invr,A)
  end function invRAR

  function RAinvR(R,A) result(m)
    real(double) :: m(3,3),R(3,3),A(3,3),tmp1(3,3),tmp3(3,3)
    integer :: i,j
    real(double),parameter :: eps=1.0d-7

    call check_orthogonal(R,eps)
    do i = 1, 3
      do j = 1, 3
        tmp1(i,j) = R(j,i)
      enddo
    enddo
    call matmatn(3,R(1,1),A(1,1),tmp3(1,1))
    call matmatn(3,tmp3(1,1),tmp1(1,1),m(1,1))
  end function RAinvR

  function get_p4_long(v1,v2,v3,r,theta,phi) result(c)
    real(double) :: v1(3),v2(3),v3(3),r,theta,phi
    real(double) :: p(3,3),newrtp(3),c(3)
    p(1:3,1) = v1(1:3)
    p(1:3,2) = v2(1:3)
    p(1:3,3) = v3(1:3)
    newrtp(1:3) = (/r,theta,phi/)
    c = get_p4(p(1,1),newrtp(1))
  end function get_p4_long

  function get_p4(p,newrtp) result(c)
    real(double) :: newrtp(3),transO(3,3),t1(3),t2(3),c(3),Ry(3,3),Rz(3,3),&
                    rtp(3),v(3),o1(3),o2(3,3),o3(3,3),p(3,3),q(3,3),r(3,3)
    integer :: i,j

    o1(1:3) = p(1:3,3)

    do i = 1, 3
      q(1:3,i) = p(1:3,i) - o1(1:3)
    enddo

    v(1:3) = q(1:3,2)
    rtp = rthetaphi(v(1))
    Ry = rotmat_y(-rtp(2)*rad2deg)
    Rz = rotmat_z(-rtp(3)*rad2deg)
    call matmatn(3,Ry(1,1),Rz(1,1),O2(1,1))

    call matmatn(3,O2(1,1),q(1,1),r(1,1))

    v(1:3) = r(1:3,1)

    rtp = rthetaphi(v(1))
    O3 = rotmat_z(-rtp(3)*rad2deg)

    v = sppolar2cart(newrtp)

    do i = 1, 3
      do j = 1, 3
        TransO(i,j) = O3(j,i)
      enddo
    enddo
    call matvecn(3,TransO(1,1),v(1),t1(1))
    do i = 1, 3
      do j = 1, 3
        TransO(i,j) = O2(j,i)
      enddo
    enddo
    call matvecn(3,TransO(1,1),t1(1),t2(1))
    c(1:3) = t2(1:3) + O1(1:3)
  end function get_p4

  subroutine get_bondlength_angle_torsion_angle(p,res)
    real(double) :: O2(3,3),Ry(3,3),Rz(3,3),rtp(3),v(3),p(3,4),q(3,4),r(3,4),res(3),s(3,4),O3(3,3)
    integer :: i

    do i = 1, 4
      q(1:3,i) = p(1:3,i) - p(1:3,3)
    enddo

    v(1:3) = q(1:3,2)
    rtp = rthetaphi(v(1))
    Ry = rotmat_y(-rtp(2)*rad2deg)
    Rz = rotmat_z(-rtp(3)*rad2deg)
    call matmatn(3,Ry(1,1),Rz(1,1),O2(1,1))

    do i = 1, 4
      call matvecn(3,O2(1,1),q(1,i),r(1,i))
    enddo

    do i = 1,4

    enddo

    v(1:3) = r(1:3,1)
    rtp = rthetaphi(v(1))

    O3 = rotmat_z(-rtp(3)*rad2deg)
    do i = 1, 4
      call matvecn(3,O3(1,1),r(1,i),s(1,i))
    enddo

    do i = 1,4

    enddo

    v(1:3) = s(1:3,4)
    rtp = rthetaphi(v(1))
    res(1:3) = (/rtp(1),rtp(2)*rad2deg,rtp(3)*rad2deg/)
  end subroutine get_bondlength_angle_torsion_angle

  function getvp(qq,v) result (c)
    real(double) :: f(3),qq(3,3),v(3),r(3),mat(3,3),Ry(3,3),&
      Rz(3,3),q(3,4),c(3),p(3)
    integer :: i,j
    do i = 1, 3
      q(1:3,i) = qq(1:3,i)
    enddo
    q(1:3,4) = v(1:3)

    f(1:3) = q(1:3,1)
    do i = 1, 4
      do j = 1, 3
        q(j,i) = q(j,i) - f(j)
      enddo
    enddo

    r = q(1:3,2)
    p = rthetaphi(r)

    Rz = rotmat_z(-p(3)*rad2deg)
    Ry = rotmat_y(-p(2)*rad2deg)
    call matmatn(3,Ry,Rz,mat)
    do i = 1, 4
      call matvecn(3,mat,q(1,i),r)
      q(1:3,i) = r
    enddo

    r = q(1:3,3)
    p = rthetaphi(r)
    Rz = rotmat_z((pi/2.0-p(3))*rad2deg)
    do i = 1, 4
      call matvecn(3,Rz,q(1,i),r)
      q(1:3,i) = r
    enddo

    Ry = rotmat_y(0.5d0*pi*rad2deg)
    do i = 1, 4
      call matvecn(3,Ry,q(1,i),r)
      q(1:3,i) = r
    enddo
    c = q(1:3,4)
  end function getvp

  subroutine unitvec(dir)
    real(double) :: dir(3),mag
    mag = vecmag3(dir(1))
    if(abs(mag) < 1.0D-10) then
      write(*,'(A)') 'vector length too small.'
    endif
    dir(1:3) = dir(1:3)/mag
  end subroutine unitvec

  function random()
    real(double)::random
    call random_number(random)
  end function random

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

  subroutine sub_min_image_dist(frac1,frac2,a,p1,p2,min_image_dist)
    real(double) :: min_image_dist,frac1(3),frac2(3),a(3,3),y(3,2)
    real(double) :: v1(3),v2(3),mindist,newdist,newfrac(1:3)
    integer :: i,j,k
    real(double) :: p1(3),p2(3)

    y(1:3,1) = frac1(1:3)
    y(1:3,2) = frac2(1:3)

    call foldtofirstzone3(y(1:3,1))
    call foldtofirstzone3(y(1:3,2))
    mindist = 1.0d10
    call frac2abs(a,y(1:3,1),v1)

    p1 = v1
    do i = -1, 1
      do j = -1, 1
        do k = -1, 1
          newfrac(1:3) = (/i,j,k/) + y(1:3,2)
          call frac2abs(a,newfrac,v2)

          newdist = vecmag3(v1-v2)
          if(newdist < mindist) then
            mindist = newdist
            p2 = v2
          endif
        enddo
      enddo
    enddo
    min_image_dist = mindist
  end subroutine sub_min_image_dist

  subroutine getnearind(frac1,frac2,a,mini,minj,mink)
    real(double) :: frac1(3),frac2(3),a(3,3),y(3,2)
    real(double) :: v1(3),v2(3),mindist,newdist,newfrac(1:3)
    integer :: i,j,k,mini,minj,mink

    y(1:3,1) = frac1(1:3)
    y(1:3,2) = frac2(1:3)
    mindist = 1.0d10

    call frac2abs(a,y(1:3,1),v1)

    do i = -1, 1
      do j = -1, 1
        do k = -1, 1

          newfrac(1:3) = (/i,j,k/) + y(1:3,2)
          call frac2abs(a,newfrac,v2)

          newdist = vecmag3(v1-v2)
          if(newdist < mindist) then
            mindist = newdist

            mini = i
            minj = j
            mink = k
          endif
        enddo
      enddo
    enddo
  end subroutine getnearind

  subroutine mole_anim_xcrys(xcrysdenfile,t,conf)
    integer :: conf
    type(supercell) :: t(:)
    character(len=*) :: xcrysdenfile

    integer :: xcrysdenu = 50,i,j
    open(unit=xcrysdenu,file=trim(xcrysdenfile),status='replace')
    write(xcrysdenu,'(A,I7)') 'ANIMSTEPS',conf
    do i = 1, conf
      write(xcrysdenu,'(A,I7)') 'ATOMS',i
      do j = 1, t(i)%n
        write(xcrysdenu,'(I4,3F20.10)') t(i)%at(j)%z,t(i)%at(j)%ac(1:3)
      enddo
    enddo
    close(xcrysdenu)
  end subroutine mole_anim_xcrys

  subroutine periodic_anim_xcrys(xfu,xname,sarray,nstr)
    character(len=*) :: xname
    integer :: nstr
    type(supercell) :: sarray(nstr)
    integer :: xfu,i,j
    real(double) :: x(3)

    open(unit=xfu,file=trim(xname),status='replace')
    write(xfu,'(A,I4)') 'ANIMSTEPS ',nstr
    write(xfu,'(A)') 'CRYSTAL'

    do j = 1, nstr
      write(xfu,'(A)') 'PRIMVEC'
      do i = 1, 3
        write(xfu,'(3f20.10)') sarray(j)%a(1:3,i)
      enddo
      write(xfu,'(A)') 'CONVVEC'
      do i = 1, 3
        write(xfu,'(3f20.10)') sarray(j)%a(1:3,i)
      enddo
      write(xfu,'(A)') 'PRIMCOORD  '//trim(N2str(j))
      if(sarray(j)%n < 1) then
        write(*,*) 'sarray(',j,')%n = ',sarray(j)%n
        write(*,'(A)') 'err: invalid n in supercell_2_xcrysden'
        stop 1
      endif
      write(xfu,'(2I4)') sarray(j)%n, 1
      do i = 1, sarray(j)%n
        call frac2abs(sarray(j)%a,sarray(j)%at(i)%f(1:3),x(1:3))
        write(xfu,'(I4,3F20.10)') sarray(j)%at(i)%z,x(1:3)
      enddo
    enddo
    close(xfu)
  end subroutine periodic_anim_xcrys

  subroutine molecule_xcrys_put(xcrysdenu,xcrysdenfile,s)
    type(supercell) :: s
    character(len=*) :: xcrysdenfile
    integer :: xcrysdenu
    integer :: i

    open(unit=xcrysdenu,file=trim(xcrysdenfile),status='replace')
    write(xcrysdenu,'(A)') 'ATOMS'
    do i = 1, s%n
      write(xcrysdenu,'(I4,3F20.10)') s%at(i)%z,s%at(i)%ac(1:3)
    enddo
    close(xcrysdenu)
  end subroutine molecule_xcrys_put

  subroutine molecule_xyz_put(xyzu,xyzfile,s)
    integer :: xyzu
    type(supercell) :: s
    character(len=*) :: xyzfile
    integer :: i

    open(unit=xyzu,file=trim(xyzfile),status='replace')
    write(xyzu,'(I6)') s%n
    write(xyzu,'(A)') 'comment'
    do i = 1, s%n
      write(xyzu,'(A,A,3F20.10)') CatS(s%at(i)%z),' ',s%at(i)%ac(1:3)
    enddo
    close(xyzu)

  end subroutine molecule_xyz_put

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

    call GetLenAng(s%a,len6(1))
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

  subroutine supercell_2_cont_vasp(originfile,contfile,s)
    type(supercell) :: s
    character(len=*) :: contfile,originfile
    integer,parameter :: contfileu=50
    integer,parameter :: poscaru=51
    integer :: i
    character(len=DSL) :: Line1,Line6

    open(unit=poscaru,file=trim(originfile),status='old',action='read')
    read(poscaru,DSF) Line1
    read(poscaru,DSF) Line6
    read(poscaru,DSF) Line6
    read(poscaru,DSF) Line6
    read(poscaru,DSF) Line6
    read(poscaru,DSF) Line6
    close(poscaru)

    open(unit=contfileu,file=trim(contfile),status='replace')

    write(contfileu,DSF) Line1
    write(contfileu,'(F6.3)') 1.0d0
    write(contfileu,'(3F18.12)') s%a(1:3,1)
    write(contfileu,'(3F18.12)') s%a(1:3,2)
    write(contfileu,'(3F18.12)') s%a(1:3,3)

    write(contfileu,DSF) Line6
    write(contfileu,'(A)') 'Direct'
    do i = 1, s%n
      write(contfileu,'(3F18.12)') s%at(i)%f(1:3)
    enddo
    close(contfileu)
  end subroutine supercell_2_cont_vasp

  function atmname_to_z(atmname)
    integer :: i,atmname_to_z
    logical :: found
    character(len=*) :: atmname

    found = .false.
    zloop: do i = 1, zordmax
      if(trim(atmname)==trim(ats(i))) then
        found = .true.
        atmname_to_z = i
        exit zloop
      endif
    enddo zloop
    if(.not. found) then
      write(*,*) 'atmname is ',trim(atmname)
      write(*,*) 'atom name not found.'
      stop 1
    endif
  end function atmname_to_z

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

  subroutine fields(str,delimit,strarr)
    character(len=*) :: str
    character :: delimit
    character(len=DSL) :: lstr
    integer :: ind
    logical :: done
    integer :: n
    type(strarrtype) :: strarr
    integer :: i

    call onespace(DSL,str,lstr)

    done = .false.

    n = 0
    strarr%n = 0
    do
      if(done) exit

      if(len(trim(lstr)) == 0) then
        done = .true.
      else
        ind = scan(lstr,delimit)

        if(ind < 1) then

          write(*,*) 'warning: problem.'
        else
          n = n+1
          strarr%s(n) = lstr(1:ind)
          strarr%n = n
        endif
        lstr = lstr(ind+1:DSL)
      endif
    enddo
    if(strarr%n == 0) then
      write(*,'(A)') 'err: no field.'
      stop 1
    else
      do i = 1, n

      enddo
    endif
  end subroutine fields

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
    call GetLenAng(sc%a,sc%la(1))
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

  function new_permute(s,n,m,nc)
    integer :: n,m,nc,i,c1,s(n),ind,new_permute

    new_permute = 1
    if(nc == 0) then
      do i = n, n+1-m, -1
        s(i) = 1
      enddo
      nc = nc + 1
    else

      ind = 0
      do i = 2, n
        if(s(i) == 1 .and. s(i-1) == -1) then
          ind = i
        endif
      enddo
      if(ind == 0) then
        new_permute = 0
      else

        c1 = 0
        do i = ind+1, n
          if(s(i) == 1) then
            c1 = c1 + 1
          endif
        enddo
        s(ind-1) = 1
        s(ind) = -1
        do i = ind+1,n-c1
          s(i) = -1
        enddo
        do i = n+1-c1,n
          s(i) = 1
        enddo
        nc = nc + 1
      endif
    endif
  end function new_permute

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

  function cohb(z)
    integer :: z
    real(double) :: cohb
    if(z < 0 .or. z .ge. zordmax) then
      write(*,*) 'z = ',z
      write(*,'(A)') 'z is too large'
    endif
    cohb = COHBOFATOM(z)

    if(abs(cohb) > 1d80) then
      write(*,*) 'z,cohb=',z,cohb
      write(*,'(A)') 'err: coherent b is not assigned yet.'
      stop 1
    endif
  end function cohb

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

  subroutine fill_in(s1,s2,nscan,orientationchoice)
    integer :: nscan,i,j,k,q,ind
    type(supercell) :: s1,s2,temp
    real(double) :: L(6),f(3),ac(3),v(3),p(3)
    real(double),parameter :: epsi=1d-4
    integer,parameter :: StandardOrientation=1
    integer,parameter :: PreserveOrientation=2
    integer :: orientationchoice

    s2%n = 0

    temp%a(1:3,1:3) = s2%a(1:3,1:3)
    call real2recip(temp%a,temp%b)

    if(orientationchoice == StandardOrientation) then

      write(*,'(A)') 'We rotate the supercell so that it has a standard orientation, i.e., a1 is along x, a2 is on the xy plane, and a3 is free'
      call GetLenAng(temp%a,L)
      call getreallatt(L,s2%a)
      call real2recip(s2%a,s2%b)
    else if(orientationchoice == PreserveOrientation) then
      write(*,*) 'We preserve the supercell orientation.'
      call real2recip(s2%a,s2%b)
    else
      write(*,*) 'orientationchoice = ',orientationchoice
      write(*,'(A)') 'err: unrecognized orientation choice.'
      stop 1
    endif

    do i = -nscan, nscan
      do j = -nscan, nscan
        do k = -nscan, nscan
          p = i*s1%a(1:3,1) + j*s1%a(1:3,2) + k*s1%a(1:3,3)
          do q = 1, s1%n
            call frac2abs(s1%a,s1%at(q)%f,ac(1:3))
            v = p + ac

            call abs2frac(temp%b,v,f)

            if(f(1) > 0d0-epsi .and. f(1) < 1d0-epsi .and. f(2) > 0d0-epsi .and. f(2) < 1d0-epsi &
              .and. f(3) > 0d0-epsi .and. f(3) < 1d0-epsi) then

              write(*,*) 'accepted frac coord: ',f(1:3)
              s2%n = s2%n + 1
              ind = s2%n
              s2%at(ind)%z = s1%at(q)%z
              call frac2abs(s2%a,f,s2%at(ind)%ac)
              call abs2frac(s2%b,s2%at(ind)%ac,s2%at(ind)%f)
            endif
          enddo
        enddo
      enddo
    enddo

    if(abs(s2%n - s1%n*det3(s2%a)/det3(s1%a)) > 1d-8) then
      write(*,*) 'Something is wrong.'
      write(*,*) 's1%n,s2%n=',s1%n,s2%n
      write(*,*) 'det3(s1%a) = ',det3(s1%a)
      write(*,*) 'det3(s2%a) = ',det3(s2%a)
      write(*,*) 'volume ratio = det3(s2%a)/det3(s1%a) = ',det3(s2%a)/det3(s1%a)
      write(*,'(A)') 'err: Simple scaling fails.'
      stop 1
    else
      write(*,*) 'fill_in: s1%n,s2%n=',s1%n,s2%n
    endif
  end subroutine fill_in

  subroutine tryput(z,x,s,newatom)
    integer :: z,i,j,k,m
    real(double) :: dis,x(3),y(3),r1(3),r2(3),f(3),q(3)
    type(supercell) :: s
    logical :: found,newatom

    call abs2frac(s%b,x,y)
    do i = 1, 3
      call modulo1(y(i))
    enddo

    call frac2abs(s%a,y,r1)
    found = .false.
    do i = 1, s%n
      if(found) exit
      f = s%at(i)%f

      do j = -1, 1
        do k = -1, 1
          do m = -1, 1
            q = (/j,k,m/) + f
            call frac2abs(s%a,q,r2)
            dis = vecmag3(r1-r2)
            if(dis < dis_tol_tryput) then
              found = .true.
              if(s%at(i)%z /= z) then
                write(*,*) 'i,s%at(i)%z=',i,s%at(i)%z
                write(*,*) 'z = ',z
                write(*,'(A)') 'err: problem. The identity of atoms are not the same.'
                stop 1
              endif
            endif
          enddo
        enddo
      enddo
    enddo
    if(.not. found) then
      s%n = s%n+1
      s%at(s%n)%f = y
      s%at(s%n)%z = z
    endif
    newatom = .not. found
  end subroutine tryput

  subroutine smallrigid_trans(s)
    type(supercell) :: s
    integer :: i,k
    real(double),parameter :: tol=1d-8
    real(double) :: x,maxmag,delta,absd

    do i = 1, 3
      maxmag = -1d10
      do k = 1, s%n
        x = s%at(k)%f(i)
        if(x < 0d0 .or. x >= 1d0) then
          write(*,'(A)') 'err: not  0 <= x < 1.'
          stop 1
        endif
        delta = x - 1d0
        absd = abs(delta)

        if(absd <= tol .and. absd > maxmag) then

          maxmag = absd
        endif
      enddo
      if(maxmag > 0) then

        do k = 1, s%n
          s%at(k)%f(i) = s%at(k)%f(i) + maxmag
          if(s%at(k)%f(i) >= 1.0d0) then
            s%at(k)%f(i) = s%at(k)%f(i) - 1.0d0
          endif
          if(abs(s%at(k)%f(i)-1.0) < 1d-12) then
            write(*,'(A)') 'err: rare cases happen in smallrigid_trans.'
            stop 1
          endif
        enddo
      endif
    enddo
  end subroutine smallrigid_trans

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

  subroutine full_jacob(rtp,B)
    real(double) :: rtp(3),B(3,3)
    real(double) :: t,p,r,ct,st,cp,sp
    r = rtp(1)
    t = rtp(2)
    p = rtp(3)
    ct = cos(t)
    st = sin(t)
    cp = cos(p)
    sp = sin(p)
    b(1,1) = st*cp
    b(2,1) = ct*cp/r
    b(3,1) = -r*sp/(r*r*st)
    b(1,2) = st*sp
    b(2,2) = ct*sp/r
    b(3,2) = r*cp/(r*r*st)
    b(1,3) = ct
    b(2,3) = -st/r
    b(3,3) = zero
  end subroutine full_jacob

  subroutine confirm_ortho(s)
    type(supercell) :: s

    real(double) :: devia

    call GetLenAng(s%a,s%la)
    devia = abs((s%la(4)-90d0)**2 + (s%la(5)-90d0)**2)
    if(devia > 1d-10) then
      write(*,*) 'devia = ',devia
      write(*,*) 's%la(4:6)=',s%la(4:6)
      write(*,'(A)') 'err: location 1 : s%la(4:6) must be 90 degs'
      stop 1
    endif
  end subroutine confirm_ortho

  subroutine merge_boxes(box1,box2,mbox)
    integer :: i,n
    type(supercell) :: box1,box2,mbox
    real(double) :: x(3)
    real(double) :: devia

    call GetLenAng(box1%a,box1%la)

    devia = abs((box1%la(4)-90d0)**2 + (box1%la(5)-90d0)**2)
    if(devia > 1d-10) then
      write(*,*) 'devia = ',devia
      write(*,*) 'box1%la(4:6)=',box1%la(4:6)
      write(*,'(A)') 'err: location 2: box1%la(4:6) must be 90 degs'
      stop 1
    endif
    call GetLenAng(box2%a,box2%la)

    devia = abs((box2%la(4)-90d0)**2 + (box2%la(5)-90d0)**2)
    if(devia > 1d-10) then
      write(*,*) 'devia = ',devia
      write(*,*) 'box2%la(4:6)=',box2%la(4:6)
      write(*,'(A)') 'err: location 3: box2%la(4:6) must be 90 degs'
      stop 1
    endif

    devia = abs((box2%la(6)-box1%la(6))**2d0)
    if(devia > 1d-10) then
      write(*,*) 'devia = ',devia
      write(*,*) 'box1%la(6)=',box1%la(6)
      write(*,*) 'box2%la(6)=',box2%la(6)
      write(*,'(A)') 'err: location 4: they must be the same.'
      stop 1
    endif

    if(abs(box1%la(1)-box2%la(1)) > 1d-10) then
       write(*,*) 'box1%la(1)=',box1%la(1)
       write(*,*) 'box2%la(1)=',box2%la(1)
       write(*,'(A)') 'err: they must be equal.'
       stop 1
    endif
    if(abs(box1%la(2)-box2%la(2)) > 1d-10) then
       write(*,*) 'box1%la(2)=',box1%la(2)
       write(*,*) 'box2%la(2)=',box2%la(2)
       write(*,'(A)') 'err: they must be equal.'
       stop 1
    endif
    mbox%la(1) = box1%la(1)
    mbox%la(2) = box1%la(2)
    mbox%la(3) = box1%la(3) + box2%la(3)

    mbox%la(4:6) = box1%la(4:6)
    call getreallatt(mbox%la,mbox%a)
    write(*,*) 'mbox%a(1:3,1) = ',mbox%a(1:3,1)
    write(*,*) 'mbox%a(1:3,2) = ',mbox%a(1:3,2)
    write(*,*) 'mbox%a(1:3,3) = ',mbox%a(1:3,3)
    call real2recip(mbox%a,mbox%b)

    n = box1%n
    do i = 1, n
      call frac2abs(box1%a,box1%at(i)%f,x)
      call abs2frac(mbox%b,x,mbox%at(i)%f)
      mbox%at(i)%z = box1%at(i)%z
    enddo
    do i = 1, box2%n
      call frac2abs(box2%a,box2%at(i)%f,x)
      x(3) = x(3) + box1%la(3)
      call abs2frac(mbox%b,x,mbox%at(n+i)%f)
      mbox%at(n+i)%z = box2%at(i)%z
    enddo
    mbox%n = n+box2%n
  end subroutine merge_boxes

  subroutine CopyBox(b1,b2)
    type(supercell) :: b1,b2
    integer :: i

    b2%commentline = b1%commentline
    b2%fracabs = b1%fracabs
    b2%la(1:6) = b1%la(1:6)
    b2%a(1:3,1:3) = b1%a(1:3,1:3)
    b2%b(1:3,1:3) = b1%b(1:3,1:3)
    b2%vol = b1%vol
    b2%n = b1%n
    b2%nsp = b1%nsp
    do i = 1, b1%n
      b2%at(i)%z = b1%at(i)%z
      b2%at(i)%f(1:3) = b1%at(i)%f(1:3)
      b2%at(i)%ac(1:3) = b1%at(i)%ac(1:3)
    enddo
  end subroutine CopyBox

  function gglrand()
    real(double) :: gglrand
    real(double),parameter :: m=2147483647d0,a=16807d0
    real(double),save :: r=10d0
    real(double) :: rr
    integer :: reminder

    rr = (a*r+ 0d0)
    reminder = nint((rr/m - int(rr/m))*m)
    r = reminder

    gglrand = reminder/m
  end function gglrand

  subroutine absdir2intarr(absdir,globalRNG)
    character(len=*) :: absdir
    type(RNGspec) :: globalRNG
    integer :: L,i

    L = len(trim(absdir))

    globalRNG%WDLEN = L
    do i = 1, L
      globalRNG%WDIntArr(i) = ichar(absdir(i:i))

    enddo
  end subroutine absdir2intarr

  subroutine read_pg_charac_table(qtablepath,tu,table,pg)
    character(len=*) :: qtablepath
    integer :: tu
    character(len=*) :: table
    type(onepointgroup) :: pg
    logical :: found
    character(len=DSL) :: tmp1,tmp2,keyword
    integer :: i,s,nc
    character(len=100) :: str1,str2

    found=.false.
    write(*,*)
    write(*,*) 'Point group: '//trim(pg%pgname)
    write(*,'(A)') 'To open '//trim(qtablepath)//'/'//trim(table)//' for reading...'
    open(unit=tu,file=trim(qtablepath)//'/'//trim(table),status='old',action='read')
    read(tu,DSF) tmp1

    do
      if(found) exit
      read(tu,DSF) tmp1
      s = scan(tmp1,' ')
      keyword = tmp1(1:s)

      if(trim(keyword)=='END') then
        write(*,'(A)') 'I cannot find the point group '//trim(pg%pgname)
        write(*,'(A)') 'err: I cannot find point group'
        stop 1
      endif
      tmp2 = tmp1(s+1:DSL)
      s= scan(tmp2,' ')
      tmp2 = tmp2(1:s)

      if(trim(tmp2) == trim(pg%pgname)) then
        found = .true.
      else

      endif
    enddo
    read(tu,*) tmp2,pg%nop,pg%nc
    if(trim(tmp2) /= 'NumOperAndNumClasses') then
      write(*,*) '2: Keyword NumOperAndNumClasses not found.'
      write(*,'(A)') 'err: Check pt-group-character-tables.'
      stop 1
    endif
    write(*,'(A,2I10)') trim(tmp2),pg%nop,pg%nc
    nc = pg%nc
    read(tu,*) tmp2
    if(trim(tmp2) /= 'ClassSize') then
      write(*,*) 'Keyword ClassSize not found.'
      write(*,'(A)') 'err: Check pt-group-character-tables.'
      stop 1
    endif
    read(tu,*) pg%classsize(1:nc)
    write(*,'(A)') 'Class sizes:'
    write(*,'(100I10)') pg%classsize(1:nc)

    read(tu,'(A200)') tmp2
    write(*,'(A)') 'tmp is '
    write(*,'(A200)') tmp2
    pg%nc = nc
    do i = 1, nc
      read(tu,*) pg%charac_mat(i,1:nc),str1,str2
      if(trim(str1) /= '||') then
        write(*,*) 'fix the character table.'
        write(*,'(A)') 'err: fix the character table.'
        stop 1
      endif
      pg%irrep(i) = trim(str2)
      write(*,*) 'irrep is '//trim(pg%irrep(i))
      write(*,'(100F10.3)') pg%charac_mat(i,1:nc)
    enddo
    close(unit=tu)
    write(*,'(A)') 'File '//trim(table)//' is closed.'
  end subroutine read_pg_charac_table

  subroutine read_pg_operations(qtablepath,pgu,sg)
    character(len=*) :: qtablepath
    type(onesg) :: sg
    integer :: pgu
    logical :: found
    character(len=DSL) :: tmp1,tmp2
    integer :: s,nop
    character(len=DSL) :: pgkeyword,numoperkeyword,pglabel
    integer :: i,j,iconfirm

    open(unit=pgu,file=trim(qtablepath)//'/table-pg',status='old',action='read')

    read(pgu,*)
    found = .false.
    do
      if(found) exit
      read(pgu,DSF) tmp1
      if(trim(tmp1) /= '') then
        s = scan(tmp1,' ')
        pgkeyword = tmp1(1:s-1)
        tmp2 = tmp1(s+1:DSL)
        s = scan(tmp2,' ')
        pglabel = tmp2(1:s-1)

      endif
      if(trim(pgkeyword)=='END') then
        write(*,*) 'PG under consideration is '//trim(sg%pg)
        write(*,'(A)') 'err: cannot find the PG label from the table'
        stop 1
      endif
      if(trim(pglabel)==trim(sg%pg)) then
        found=.true.
      else

      endif
    enddo
    if(.not. found) then
      write(*,*) 'trim(sg%pg) is '//trim(sg%pg)
      write(*,'(A)') 'err: not found.'
      stop 1
    else

    endif

    read(pgu,DSF) tmp1
    read(tmp1,*) numoperkeyword, nop
    if(trim(numoperkeyword) /= 'NumOper') then
      write(*,*) 'numoperkeyword is '//trim(numoperkeyword)
      write(*,'(A)') 'err: NumOper must be found.'
      stop 1
    endif
    write(*,'(A)') 'Point group '//trim(sg%pg)//' has '//trim(N2str(nop))//' operations'
    sg%PGNumOper = nop
    do i = 1, nop
      read(pgu,*) iconfirm
      if(iconfirm/= i) then
        write(*,*) 'iconfirm,i=',iconfirm,i
        write(*,'(A)') 'err: iconfirm is not equal to i'
        stop 1
      else

      endif
      sg%op(1,1:4,i) = (/one,zero,zero,zero/)
      do j = 2, 4
        read(pgu,*) sg%op(j,2:4,i)
      enddo
    enddo
    close(pgu)
  end subroutine read_pg_operations

  function pg_id_map(d,t)
    integer :: pg_id_map,d,t

    pg_id_map = d*pg_uuid_L4 + abs(t)*pg_uuid_L3 + t*pg_uuid_L2 + t*d
  end function pg_id_map

  subroutine PG_NAME_init()
    PG_NAME(1) =  trim("C1")
    PG_NAME(2) =  trim("Ci")
    PG_NAME(3) =  trim("C2")
    PG_NAME(4) =  trim("Cs")
    PG_NAME(5) =  trim("C2h")
    PG_NAME(6) =  trim("D2")
    PG_NAME(7) =  trim("C2v")
    PG_NAME(8) =  trim("D2h")
    PG_NAME(9) =  trim("C4")
    PG_NAME(10) = trim("S4")
    PG_NAME(11) = trim("C4h")
    PG_NAME(12) = trim("D4")
    PG_NAME(13) = trim("C4v")
    PG_NAME(14) = trim("D2d")
    PG_NAME(15) = trim("D4h")
    PG_NAME(16) = trim("C3")
    PG_NAME(17) = trim("S6")
    PG_NAME(18) = trim("D3")
    PG_NAME(19) = trim("C3v")
    PG_NAME(20) = trim("D3d")
    PG_NAME(21) = trim("C6")
    PG_NAME(22) = trim("C3h")
    PG_NAME(23) = trim("C6h")
    PG_NAME(24) = trim("D6")
    PG_NAME(25) = trim("C6v")
    PG_NAME(26) = trim("D3h-type1")
    PG_NAME(27) = trim("D6h")
    PG_NAME(28) = trim("T")
    PG_NAME(29) = trim("Th")
    PG_NAME(30) = trim("O")
    PG_NAME(31) = trim("Td")
    PG_NAME(32) = trim("Oh")
  end subroutine PG_NAME_init

  function get_pg_index(pg_id)
    integer :: pg_id,get_pg_index
    logical :: found
    integer :: i,ind

    found = .false.
    ind = -1
    do i = 1, 32
      if(.not. found .and. pg_id == PG_UUID(i)) then
        found = .true.
        ind = i
      endif
    enddo
    if(.not. found) then
      write(*,'(A)') 'err: get_pg_index: UUID not found.'
      stop 1
    endif
    get_pg_index = ind
  end function get_pg_index

  subroutine check_SG_closure(sg)
    type(onesg) :: sg
    integer :: nop
    real(double),allocatable :: op(:,:,:),origmat3(:,:,:),trace(:),pg_det3(:)
    integer :: i,j,k
    real(double) :: tmpmat(4,4),p(3),diffmat(4,4),diffmat3(3,3),mat1(3,3),mat2(3,3),mat3(3,3)
    real(double),parameter :: eps = 1.0d-8
    integer,allocatable :: sudoku(:,:),classes(:)
    logical :: found
    real(double) :: matmag,sumj, targetsum
    real(double),allocatable :: sortv(:)
    real(double) :: RG,mattrace,detm
    integer,allocatable :: ig(:),indarr(:)
    integer :: NG,abs_op,det_int,trace_int
    integer :: pg_id,pg_index
    logical :: print_op
    real(double) :: diffmag

    write(*,*)
    write(*,*) 'Checking space group information.'
    write(*,*) 'Name of SG is '//trim(sg%sglabel)
    write(*,*) 'Name of PG is '//trim(sg%pg)
    nop = sg%nop
    write(*,*) 'Number of operations= ',trim(N2str(nop))

    allocate(op(4,4,nop))

    print_op = .false.
    if(print_op) then
      write(*,'(A)') 'In check_SG_closure (for cubic SGs, the order might have'
      write(*,'(A)') 'been ordered according to conjugacy classes):'
    endif

    do i = 1, nop
      op(1:4,1:4,i) = sg%op(1:4,1:4,i)

      if(print_op) then
        write(*,*) 'i = ', i

        diffmag = abs(op(1,1,i)-one) + abs(op(1,2,i)) + abs(op(1,3,i)) + abs(op(1,4,i))
        if(diffmag > 1.0d-16) then
          write(*,'(A)') 'err: first row must be (1,0,0,0)'
          stop 1
        endif
        write(*,'(4F8.3)') op(2,1:4,i)
        write(*,'(4F8.3)') op(3,1:4,i)
        write(*,'(4F8.3)') op(4,1:4,i)
      endif
    enddo

    allocate(sudoku(nop,nop))

    do i = 1, nop
      do j = 1, nop
        call matmatn(4,op(1,1,i),op(1,1,j),tmpmat(1,1))

        p(1:3) = tmpmat(2:4,1)
        call foldtofirstzone3(p(1))
        do k = 1, 3
          if(p(k) > one - eps) then
            write(*,*) 'p(k) = ', p(k)
            p(k) = zero
            write(*,*) 'reset p(k), p(k) = ',p(k)
            write(*,'(A)') 'err: does this really happen?'
            stop 1
          endif
        enddo

        tmpmat(2:4,1) = p(1:3)

        found = .false.
        do k = 1, nop
          diffmat(1:4,1:4) = sg%op(1:4,1:4,k) - tmpmat(1:4,1:4)
          matmag = matmagn(4,diffmat(1,1))
          if(matmag < eps) then
            found = .true.
            sudoku(i,j) = k
          endif
        enddo
        if( .not. found) then
          write(*,*) 'i, j= ',i,j
          write(*,*) 'tmpmat(1:4,1:4)= '
          write(*,'(4F8.3)') tmpmat(1,1:4)
          write(*,'(4F8.3)') tmpmat(2,1:4)
          write(*,'(4F8.3)') tmpmat(3,1:4)
          write(*,'(4F8.3)') tmpmat(4,1:4)
          write(*,*) 'closure test failed.'
          write(*,'(A)') 'err: closure fail.'
          stop 1
        endif
      enddo
    enddo
    write(*,*)

    do i = 1, nop

    enddo
    targetsum = nop*(nop+one)/two

    do i = 1, nop
      sumj = zero
      do j = 1, nop
        sumj = sumj + sudoku(i,j)
      enddo
      if(abs(sumj - targetsum) > eps) then
        write(*,*) 'failed: sumj, targetsum= ',sumj, targetsum
        write(*,'(A)') 'err: matrix multiplication table sudoku property failed.'
        stop 1
      endif
    enddo

    do i = 1, nop
      sumj = zero
      do j = 1, nop
        sumj = sumj + sudoku(j,i)
      enddo
      if(abs(sumj - targetsum) > eps) then
        write(*,*) 'sumj, targetsum= ',sumj, targetsum
        write(*,'(A)') 'err: sudoku failed.'
        stop 1
      endif
    enddo
    write(*,'(A)') 'Passed closure property test'

    if(nop /= sg%PGNumOper) then

      write(*,*) 'WARNING: We will use the exact number of operations for the PG.'
      write(*,'(A,I4,/,A,I4)') 'We change the number of operations in the SG= ',nop,' to number of operations in the PG= ',sg%PGNumOper
    endif
    nop = sg%PGNumOper

    allocate(trace(nop))
    allocate(pg_det3(nop))

    do i = 1, nop
      sumj = zero
      do j = 2, 4
        sumj = sumj + op(j,j,i)
      enddo
      trace(i) = sumj

    enddo

    allocate(classes(nop))
    allocate(origmat3(3,3,nop))
    do i = 1, nop
      origmat3(1:3,1:3,i) = op(2:4,2:4,i)
      pg_det3(i) = det3(origmat3(1,1,i))
    enddo

    classes(1:nop) = 0
    do i = 1, nop
      if(classes(i) /= 0) cycle
      mat1(1:3,1:3) = origmat3(1:3,1:3,i)

      do j = 1, nop

        mat2 = origmat3(1:3,1:3,j)

        mat3 = BAinvB(mat2(1,1),mat1(1,1))

        found = .false.
        do k = 1, nop
          diffmat3(1:3,1:3) = mat3(1:3,1:3)-origmat3(1:3,1:3,k)

          if(.not. found .and. matmagn(3,diffmat3(1,1)) < eps) then
            found = .true.
            if(classes(k) /= 0) then
              if(classes(k) /= i) then
                write(*,*) 'i,classes(k)= ',i,classes(k)
                write(*,'(A)') 'err: inconsistency.'
                stop 1
               endif
            else
              classes(k) = i
            endif
          endif
        enddo
        if( .not. found) then
          write(*,'(A)') 'err: cannot find companion in equivalent classes.'
          stop 1
        endif
      enddo
    enddo
    if(nop < 1) then
      write(*,*) 'nop = ',nop
      write(*,'(A)') 'err: going to hit allocation problem since the array size is zero.'
      stop 1
    endif
    allocate(sortv(nop))
    allocate(indarr(nop))
    do i = 1, nop
      sortv(i) = classes(i)
      indarr(i) = i
    enddo
    call dbl_sort(nop,sortv(1),indarr(1),1)

    do i = 1, nop

    enddo
    do i = 1, nop
      if(i == indarr(i)) then

      else

        write(*,'(A)') ' We expect the operations have been fully ordered.'
        write(*,'(A)') 'err: Search the keyword OPERATION_MAP in commod.f90.'
        stop 1
      endif
    enddo
    write(*,'(A)') 'Passed the natural ordering test for equivalent classes.'

    allocate(ig(nop+1))
    RG = sortv(1)
    NG = 1
    IG(1) = 1
    do k = 2, nop
      if(abs(sortv(k)-RG) > 1.0d-4) then
        RG = sortv(k)
        NG = NG + 1
        IG(NG) = k
      endif
    enddo
    ig(ng+1) = nop+1

    write(*,*)
    write(*,'(A)') 'Equivalence classes of the point group: '//trim(sg%pg)

    pg_id = ng*PG_UUID_L5
    write(*,'(A)') '--------- classes begin -----'
    do i = 1, ng
      abs_op = indarr(ig(i))
      detm = det3( origmat3(1,1, abs_op ))
      mattrace = tracen(3, origmat3(1,1,abs_op))
      write(*,'(A,I4,I4,2F8.4)') 'cluster, size of cluster, det, trace=',i,ig(i+1)-ig(i),detm,mattrace
      det_int = nint(detm)
      trace_int = nint(mattrace)
      call report_operation_type(det_int,trace_int)
      pg_id  = pg_id + pg_id_map(det_int,trace_int)
    enddo
    write(*,'(A)') '--------- classes end -----'
    pg_index = get_pg_index(pg_id)
    deallocate(op)
    deallocate(sudoku)
    deallocate(trace)
    deallocate(pg_det3)
  end subroutine check_SG_closure

  subroutine read_sg(qtablepath,sg_label,sg,pgu,tu,dtu,sgtable,dsgtable)
    type(onesg) :: sg
    character(len=*) :: sg_label,qtablepath
    character(len=*) :: sgtable,dsgtable
    integer :: pgu,tu,dtu
    character(len=DSL) :: base
    character(len=DSL) :: lpart
    character(len=DSL) :: ext
    integer :: str_index
    integer :: strlen,ind,newnc,newnop
    integer :: readmode

    strlen = len(trim(sg_label))
    ind = scan(sg_label,':')

    if(ind > 0) then
      write(*,*) 'sg_label contains :'
      lpart = sg_label(1:ind-1)
      ext = sg_label(ind+1:strlen)
      write(*,*) 'lpart is '//trim(lpart)
      write(*,*) 'ext is '//trim(ext)

      str_index = index(trim(sg_label),'basic')
      if(str_index > 0) then
        write(*,*) 'found basic in the sglabel'
        base = sg_label(1:str_index-2)
        write(*,*) 'base is '//trim(base)
        call sg_table_read(qtablepath,pgu,tu,sgtable,base,sg,readmode)
        call check_SG_closure(sg)

        newnc = 1
        newnop = sg%PGNumOper
        write(*,'(/A)') 'With : case, since the label has "-basic", we set ncentering from '//trim(N2str(sg%ncentering))//' to '//trim(N2str(newnc))//', and nop from '//trim(N2str(sg%nop))//' to '//trim(N2str(newnop))
        sg%ncentering = newnc
        sg%nop = newnop

      else
        write(*,*) 'no basic in the sglabel'

        base = trim(lpart)
        call sg_table_read(qtablepath,pgu,tu,sgtable,base,sg,readmode)
        call check_SG_closure(sg)
      endif

      call update_derivsg_table(qtablepath,dtu,dsgtable,sg_label,sg,lpart,ext)

      call check_SG_closure(sg)
    else
      str_index = index(trim(sg_label),'basic')
      if(str_index > 0) then
        write(*,*) 'found basic in the sglabel'
        base = sg_label(1:str_index-2)
        write(*,*) 'base is '//trim(base)
        call sg_table_read(qtablepath,pgu,tu,sgtable,base,sg,readmode)
        if(readmode == 1 .or. readmode == 2) then

        else
          write(*,*) 'readmode =',readmode
          write(*,'(A)') 'err: impossible.'
          stop 1
        endif

        call check_SG_closure(sg)

        write(*,'(/A)') 'Without : case, since the label has "-basic", we set ncentering to '//trim(N2str(sg%ncentering))//' and nop to '//trim(N2str(sg%nop))
        sg%ncentering = 1
        sg%nop = sg%PGNumOper

      else

        base = trim(sg_label)
        call sg_table_read(qtablepath,pgu,tu,sgtable,base,sg,readmode)
        if(readmode == 1 .or. readmode == 2) then
          call check_SG_closure(sg)
        else if(readmode == 3) then
          write(*,*) 'readmode=3: we do not check the SG_closure.'
        else
          write(*,*) 'readmode =',readmode
          write(*,'(A)') 'err: readmode problem.'
          stop 1
        endif
      endif
    endif
  end subroutine read_sg

  subroutine sg_table_read(qtablepath,pgu,tu,table,targetlabel,sg,readmode)
    integer :: pgu,tu
    character(len=*) :: table,qtablepath
    character(len=*) :: targetlabel
    type(onesg) :: sg
    logical :: found
    integer :: i,j,k,iconfirm,nop
    character(len=DSL) :: sglabel
    character(len=DSL) :: sgkeyword,numoperkeyword
    character(len=DSL) :: tmp1,tmp2,tmp3,str,key
    integer :: s,v,ncentering,nopread
    character(len=DSL) :: string(3)
    real(double) :: f,nf,df
    integer :: pru
    real(double) :: shift(3)
    real(double) :: tmpMat(4,4,48)
    character(len=DSL) :: sg_f90_entry,master_sg_entry
    integer :: readmode
    integer :: labellen

    readmode = -10000

    write(*,*)
    labellen = len(trim(targetlabel))
    if(labellen == 0) then
      write(*,*) 'Empty targetlabel.'
      write(*,'(A)') 'err: Empty targetlabel.'
      stop 1
    endif
    write(*,'(A)') 'sg_table_read: targetlabel is '//trim(targetlabel)

    sg%assigned = .FALSE.

    do i = 0, 230
      do j = 1, maxSGvariant
        if(trim(SGbase(j,i)%sglabel)==trim(targetlabel)) then
          sg%sglabel = trim( SGbase(j,i)%sglabel)
          sg%PG = trim(SGbase(j,i)%PG)
          sg%ncentering = SGbase(j,i)%ncentering
          do k = 1, sg%ncentering
            sg%centering(1:3,k) = SGbase(j,i)%centering(1:3,k)
          enddo

          call read_pg_operations(qtablepath,pgu,sg)

          if(sg%PGNumOper /= SGbase(j,i)%PGNumOper) then
            write(*,*) 'sg%PGNumOper,SGbase(j,i)%PGNumOper=',sg%PGNumOper,SGbase(j,i)%PGNumOper
            write(*,'(A)') 'err: serious inconsistent problem.'
            stop 1
          endif
          sg%symmorphic = SGbase(j,i)%symmorphic
          if(sg%symmorphic == 1) then
            do k = 1, sg%PGNumOper
              sg%nsy(1:3,k) = (/zero,zero,zero/)
              sg%op(2:4,1,k) = (/zero,zero,zero/)
            enddo
          else
            do k = 1, sg%PGNumOper
              sg%nsy(1:3,k) = SGbase(j,i)%nsy(1:3,k)
              sg%op(2:4,1,k) = SGbase(j,i)%nsy(1:3,k)
            enddo
          endif
          sg%nop = SGbase(j,i)%nop
          sg%reorder_IT_op_sequence = SGbase(j,i)%reorder_IT_op_sequence
          if(sg%reorder_IT_op_sequence == 1) then
            sg%identity_map(1:sg%PGNumOper) = SGbase(j,i)%identity_map(1:sg%PGNumOper)
            sg%operation_map(1:sg%PGNumOper) = SGbase(j,i)%operation_map(1:sg%PGNumOper)
          endif
          sg%assigned = .TRUE.
          readmode = 1
        endif
      enddo
    enddo

    if(sg%assigned) then

      nop = sg%nop
      ncentering = nop/sg%PGNumOper
      if(ncentering /= sg%ncentering) then
        write(*,*) 'ncentering,sg%ncentering = ',ncentering,sg%ncentering
        write(*,'(A)') 'err: inconsistency.'
        stop 1
      endif

      if(sg%reorder_IT_op_sequence == 1) then
        do i = 1, sg%PGNumOper
          tmpMat(1:4,1:4,i) = sg%op(1:4,1:4,i)
        enddo
        do i = 1, sg%PGNumOper
          sg%op(1:4,1:4,i) = tmpMat(1:4,1:4,   sg%operation_map(i))
        enddo
        write(*,'(A)') 'In mode 1: Operations are reordered for O,Td,Oh since IT SGs equivalence classes are not contiguous'
      endif

      do i = 2, ncentering
        do j = 1, sg%PGNumOper
          sg%op(1:4,1:4,(i-1)*sg%PGNumOper+j) = sg%op(1:4,1:4,j)

          sg%op(2:4,1,(i-1)*sg%PGNumOper+j) = sg%op(2:4,1,(i-1)*sg%PGNumOper+j) + sg%centering(1:3,i)

          shift(1:3) = sg%op(2:4,1,(i-1)*sg%PGNumOper+j)
          call foldtofirstzone3(shift(1))

          sg%op(2:4,1,(i-1)*sg%PGNumOper+j) = shift(1:3)

        enddo
      enddo

    else
      write(*,'(A)') 'Do not use a SG entry in commod.f90, instead we use '//trim(table)

      write(*,'(A)') ' -------- To open SG table file '//trim(qtablepath)//'/'//trim(table)
      open(unit=tu,file=trim(qtablepath)//'/'//trim(table),status='old',action='read')

      read(tu,*)

      found=.false.

      do
        if(found) exit

        read(tu,DSF) tmp1

        if(trim(tmp1) /= '') then

          s = scan(tmp1,' ')
          sgkeyword = tmp1(1:s-1)
          tmp2 = tmp1(s+1:DSL)
          s = scan(tmp2,' ')
          sglabel = tmp2(1:s-1)

        endif

        if(trim(sgkeyword)=='END') then
          exit
        endif

        if(trim(sglabel)==trim(targetlabel)) then
          found=.true.
        else

        endif
      enddo

      if(found) then
        write(*,'(A)') 'Found SG label= '//trim(targetlabel)
      else
        write(*,'(A)') 'SG under consideration is '//trim(targetlabel)
        write(*,'(A)') 'err: in sg_table_read, cannot find SG label.'
        stop 1
      endif

      sg%PGNumOper=-1

      sg%ncentering=-1
      sg%reorder_IT_op_sequence = 0

      do
        read(tu,DSF) tmp1

        read(tmp1,*) str

        if(trim(str) == '#optional') then
          s = scan(str,' ')
          tmp2 = tmp1(s+1:DSL)

          s = scan(tmp2,' ')
          key=tmp2(1:s)

          tmp3 = tmp2(s+1:DSL)

          if(trim(key) == 'PG') then
            readmode = 2
            write(*,'(A)') 'Use a PG label to assign the matrix operations.'
            write(*,'(A)') 'PG: to process '//trim(tmp3)
            s = scan(tmp3,' ')
            sg%pg = tmp3(1:s)
            write(*,'(A)') 'Point group '//trim(sg%pg)

          else if(trim(key) == 'ncentering') then

            read(tmp3,'(I10)') v
            write(*,*) 'ncentering = ',v
            if(v < 1 .or. v > 4) then
              write(*,*) 'v = ',v
              write(*,'(A)') 'err: invalid v.'
              stop 1
            endif
            sg%ncentering = v
            do i = 1, v
              read(tu,*) tmp1,tmp2,sg%centering(1:3,i)
              iconfirm = str2N(tmp2)
              if(i /= iconfirm) then
                write(*,*) 'centering: i,iconfirm=',i,iconfirm
                write(*,'(A)') 'err: i and iconfirm must be the same.'
                stop 1
              endif
              write(*,'(A,3F10.5,A)') 'tmp1 is '//trim(tmp1)//', tmp2 is '//trim(tmp2)//', ncentering= |', sg%centering(1:3,i),'|'
            enddo

          else if(trim(key) == 'PGNumOper') then

            call read_pg_operations(qtablepath,pgu,sg)

            read(tmp3,'(I10)') v
            write(*,*) 'PGNumOper = ',v
            if(sg%PGNumOper /= v) then
              write(*,*) 'sg%PGNumOper = ',sg%PGNumOper
              write(*,*) 'v = ',v
              write(*,'(A)') 'err: serious inconsistency problem.'
              stop 1
            endif

            read(tu,*) tmp1,tmp2 !
            write(*,'(A,A)') 'tmp1 = '//trim(tmp1),', tmp2 = '//trim(tmp2)
            if(trim(tmp2) == 'symmorphic' .or. trim(tmp2) == 'nonsymmorphic') then

            else
              write(*,'(A)') 'err: the keyword must be symmorphic and nonsymmorphic.'
              stop 1
            endif

            if(trim(tmp2) == 'symmorphic') then
              sg%symmorphic = 1

              do i = 1, v
                sg%op(2:4,1,i) = (/zero,zero,zero/)
              enddo
            else
              sg%symmorphic = 0
              write(*,*) 'v = ',v
              do i = 1, v
                read(tu,*) tmp1,tmp2,sg%op(2:4,1,i)
                iconfirm = str2N(tmp2)
                if(i /= iconfirm) then
                  write(*,*) 'fractional translation: i,iconfirm=',i,iconfirm
                  write(*,'(A)') 'err: i and iconfirm must be the same.'
                  stop 1
                endif
                if(i > 48) then
                  write(*,'(A)') 'err: bound exceeded.'
                  stop 1
                endif
                sg%nsy(1:3,i) = sg%op(2:4,1,i)
                write(*,'(A,3F10.5)') 'tmp1 is '//trim(tmp1)//', tmp2 is '//trim(tmp2)//', translation vector= ', sg%op(2:4,1,i)
              enddo
            endif

          else if(trim(key) == 'reorder_IT_op_sequence') then
            if(sg%PGNumOper < 1) then
              write(*,'(A)') 'err: bad situation.'
              stop 1
            endif
            sg%reorder_IT_op_sequence = 1
            read(tu,*) tmp1, tmp2, sg%identity_map(1:sg%PGNumOper)
            if(trim(tmp1) /= '#optional' .or. trim(tmp2) /= 'identity_map') then
              write(*,*) 'tmp1= '//trim(tmp1),', tmp2= ',trim(tmp2)
              write(*,'(A)') 'err: keyword mismatched.'
              stop 1
            endif
            read(tu,*) tmp1, tmp2, sg%operation_map(1:sg%PGNumOper)
            if(trim(tmp1) /= '#optional' .or. trim(tmp2) /= 'operation_map') then
              write(*,*) 'tmp1= '//trim(tmp1),', tmp2= ',trim(tmp2)
              write(*,'(A)') 'err: keyword mismatched.'
              stop 1
            endif
          else
            write(*,*) 'keyword of |'//trim(key)//'| will not be processed.'
            write(*,'(A)') 'err: keyword not found.'
            stop 1
          endif
        else

          exit
        endif
      enddo

      if(readmode == 2) then
        if(sg%ncentering == -1 .or. sg%PGNumOper == -1) then
          write(*,*) 'sg%ncentering = ',sg%ncentering
          write(*,*) 'sg%PGNumOper = ',sg%PGNumOper
          write(*,'(A)') 'err: Serious inconsistency input.'
          stop 1
        endif
      endif

      if(sg%reorder_IT_op_sequence == 1) then
        do i = 1, sg%PGNumOper
          tmpMat(1:4,1:4,i) = sg%op(1:4,1:4,i)
        enddo
        do i = 1, sg%PGNumOper
          sg%op(1:4,1:4,i) = tmpMat(1:4,1:4,   sg%operation_map(i))
        enddo
        write(*,'(A)') 'In mode 2: Operations are reordered for O,Td,Oh since IT SGs equivalence classes are not contiguous'
      endif

      read(tmp1,*) numoperkeyword, nop

      if(trim(numoperkeyword) /= 'NumOper') then
        write(*,*) 'numoperkeyword is '//trim(numoperkeyword)
        write(*,'(A)') 'err: NumOper must be found.'
        stop 1
      endif
      write(*,'(A)') 'space group of '//trim(targetlabel)//' has '//trim(N2str(nop))//' operations'

      sg%sglabel = trim(targetlabel)
      sg%nop = nop

      if(readmode == 2) then

        if(mod(nop, sg%PGNumOper) /= 0) then
          write(*,'(A)') 'err: PGNumOper does not divide nop.'
          stop 1
        endif
        ncentering = nop/sg%PGNumOper
        if(ncentering /= sg%ncentering) then
          write(*,*) 'ncentering,sg%ncentering = ',ncentering,sg%ncentering
          write(*,'(A)') 'err: Inconsistency.'
          stop 1
        endif

        do i = 2, ncentering
          do j = 1, sg%PGNumOper
            sg%op(1:4,1:4,(i-1)*sg%PGNumOper+j) = sg%op(1:4,1:4,j)
            sg%op(2:4,1,(i-1)*sg%PGNumOper+j) = sg%op(2:4,1,(i-1)*sg%PGNumOper+j) + sg%centering(1:3,i)

            shift(1:3) = sg%op(2:4,1,(i-1)*sg%PGNumOper+j)
            call foldtofirstzone3(shift(1))
            sg%op(2:4,1,(i-1)*sg%PGNumOper+j) = shift(1:3)
          enddo
        enddo

      else
        readmode = 3
        write(*,*) 'This will be defaulted to reading operation by operation from the table'
        nopread = nop

        do i = 1, nopread
          read(tu,*) iconfirm
          if(iconfirm/= i) then
            write(*,*) 'iconfirm,i=',iconfirm,i
            write(*,'(A)') 'err: iconfirm is not equal to i'
            stop 1
          else

          endif
          sg%op(1,1:4,i) = (/one,zero,zero,zero/)
          do j = 2, 4
            read(tu,*) sg%op(j,1:4,i)
          enddo
        enddo

        return
      endif
      close(tu)
      write(*,'(A)') ' ------- Closing SG table file '//trim(table)
    endif

    if(readmode == 1 .or. readmode == 2) then
    else
      write(*,*) 'readmode = ',readmode
      write(*,'(A)') 'err: bad readmode.'
      stop 1
    endif

    write(*,*)
    sg_f90_entry='echo-sg-f90-entry.dat'

    pru = newunit()
    open(unit=pru,file=trim(sg_f90_entry),status='replace')

    write(pru,'(A)') '--------------------------------------'

    sglabel = trim(sg%sglabel)
    i = scan(sglabel,'-')
    write(pru,'(A)') '    SG = '//sglabel(1:i-1)
    write(pru,'(A)') '    Y = 1'
    write(pru,'(A)') '    SGbase(Y,SG)%sglabel    ='''//trim(sg%sglabel)//''''
    write(pru,'(A)') '    SGbase(Y,SG)%PG         ='''//trim(sg%pg)//''''
    write(pru,'(A)') '    SGbase(Y,SG)%ncentering ='//trim(N2str(sg%ncentering))
    do i = 1, sg%ncentering
      do j = 1, 3
        f = sg%centering(j,i)

        call frac2nd(f,nf,df)
        if(abs(df-one) < 1.0d-6) then

          string(j)=trim(N2str(nint(nf)))//'.0d0'
        else
          string(j)=trim(N2str(nint(nf)))//'.0d0/'//trim(N2str(nint(df)))//'.0d0'
        endif
      enddo
      write(pru,'(A,A16,A,A16,A,A16,A)') '    SGbase(Y,SG)%centering(1:3,'//trim(N2str(i))//') = (/',trim(string(1)),',',trim(string(2)),',',trim(string(3)),' /)'
    enddo
    write(pru,'(A)') '    SGbase(Y,SG)%PGNumOper  ='//trim(N2str(sg%PGNumOper))
    write(pru,'(A)') '    SGbase(Y,SG)%symmorphic ='//trim(N2str(sg%symmorphic))
    if(sg%symmorphic == 0) then
      do i = 1, sg%PGNumOper
        do j = 1, 3
          f = sg%nsy(j,i)

          call frac2nd(f,nf,df)
          if(abs(df-one) < 1.0d-6) then

            string(j)=trim(N2str(nint(nf)))//'.0d0'
          else
            string(j)=trim(N2str(nint(nf)))//'.0d0/'//trim(N2str(nint(df)))//'.0d0'
          endif
        enddo
        write(pru,'(A,A16,A,A16,A,A16,A)') '    SGbase(Y,SG)%nsy(1:3,'//trim(N2str(i))//') = (/',trim(string(1)),',',trim(string(2)),',',trim(string(3)),' /)'
      enddo
    endif

    if(sg%reorder_IT_op_sequence == 1) then
      write(pru,'(A)') '    SGbase(Y,SG)%reorder_IT_op_sequence = 1'
      tmp1=trim(N2str(sg%identity_map(1)))
      do i = 2, sg%PGNumOper
        tmp1=trim(tmp1)//','//trim(N2str(sg%identity_map(i)))
      enddo
      write(pru,'(A)') '    SGbase(Y,SG)%identity_map(1:'//trim(N2str(sg%PGNumOper))//') = (/'//trim(tmp1)//'/)'
      tmp1=trim(N2str(sg%operation_map(1)))
      do i = 2, sg%PGNumOper
        tmp1=trim(tmp1)//','//trim(N2str(sg%operation_map(i)))
      enddo
      write(pru,'(A)') '    SGbase(Y,SG)%operation_map(1:'//trim(N2str(sg%PGNumOper))//') = (/'//trim(tmp1)//'/)'
    else

    endif

    write(pru,'(A)') '    SGbase(Y,SG)%nop        ='//trim(N2str(sg%nop))
    write(pru,'(A)') '    SGbase(Y,SG)%assigned   =.TRUE.'
    write(pru,*)
    close(pru)

    master_sg_entry='echo-master-sg-entry.dat'

    pru = newunit()
    open(unit=pru,file=trim(master_sg_entry),status='replace')
    write(pru,'(A)') 'SG '//trim(sg%sglabel)
    write(pru,'(A)') '#optional PG '//trim(sg%pg)
    write(pru,'(A)') '#optional ncentering '//trim(N2str(sg%ncentering))
    do i = 1, sg%ncentering
      write(pru,'(A,3F20.15)') '#optional '//trim(N2str(i)),sg%centering(1:3,i)
    enddo
    write(pru,'(A)') '#optional PGNumOper '//trim(N2str(sg%PGNumOper))
    if(sg%symmorphic == 1) then
      write(pru,'(A)') '#optional symmorphic'
    else if(sg%symmorphic == 0) then
      write(pru,'(A)') '#optional nonsymmorphic'
      do i = 1, sg%PGNumOper
        write(pru,'(A,3F8.4)') '#optional '//trim(N2str(i)),sg%nsy(1:3,i)
      enddo
    else
      write(*,'(A)') 'err: problem in symmorphic.'
      stop 1
    endif
    if(sg%reorder_IT_op_sequence == 1) then
      write(pru,'(A)') '#optional reorder_IT_op_sequence'
      write(tmp2,'(100I3)') sg%identity_map(1:sg%PGNumOper)
      write(pru,'(A)') '#optional  identity_map '//trim(tmp2)
      write(tmp2,'(100I3)') sg%operation_map(1:sg%PGNumOper)
      write(pru,'(A)') '#optional operation_map '//trim(tmp2)
    else

    endif
    write(pru,'(A)') 'NumOper '//trim(N2str(sg%PGNumOper*sg%ncentering))

    close(pru)

  end subroutine sg_table_read

  subroutine update_derivsg_table(qtablepath,dtu,dtable,targetlabel,sg,lpart,ext)
    character(len=*) :: qtablepath
    integer :: dtu
    character(len=*) :: dtable
    character(len=*) :: targetlabel
    type(onesg) :: sg
    logical :: found
    character(len=DSL) :: line,tmp1,tmp2,sgkeyword,sglabel,ext,lpart
    integer :: s
    real(double) :: invT(3,3),blk4x4(4,4),newblk4x4(4,4),t1(3),R1(3,3),t2(3),R2(3,3)
    integer :: i,nop

    real(double) :: t3(3)
    real(double) :: origvec(3)
    integer :: j,pru
    real(double) :: f,nf,df
    character(len=DSL) :: string(3)
    character(len=DSL) :: filename
    integer :: ind,SGn
    character :: FI

    write(*,*)
    write(*,*)
    write(*,*) 'in update_derivsg_table.'

    if(len(trim(targetlabel)) == 0) then
      write(*,*) 'Empty targetlabel.'
      stop
    else
      write(*,*) 'targetlabel is '//trim(targetlabel)
    endif

    if(trim(ext) == 'stdbox' .or. trim(ext) == 'mat2' .or. trim(ext) == 'mat6') then
      write(*,*) 'lpart is '//trim(lpart)
      ind = scan(lpart,'-')
      tmp1 = lpart(1:ind-1)
      read(tmp1,'(I10)') SGn
      FI = lpart(ind+1:ind+1)
      write(*,*) 'SGn,FI = ',SGn,FI

      if(SGn >= 1 .and. SGn <= 2) then
        write(*,*) 'Triclinic case does not a chance to use stdbox'
        write(*,'(A)') 'triclinic'
      else if(SGn >= 3 .and. SGn <= 15) then
        if(trim(ext) /= 'stdbox') then
          write(*,*) 'For monoclinic, we must have stdbox case'
          write(*,'(A)') 'err: not stdbox.'
          stop 1
        endif

        s = index(lpart,'-')
        FI = lpart(s+1:s+1)
        if(FI == 'C') then
          sg%T(1:3,1:3) = MonocBaxisBaseCentredCmat(1:3,1:3)
        else
           write(*,'(A)') 'er: crash, P should call this subroutine.'
        endif
      else if(SGn >= 16 .and. SGn <= 74) then
        if(trim(ext) /= 'stdbox') then
          write(*,'(A)') 'err: not stdbox.'
          stop 1
        endif
        FI = lpart(4:4)
        if(FI == 'C') then
          sg%T(1:3,1:3) = OBaseCmat1(1:3,1:3)
        else if (FI == 'F') then
          sg%T(1:3,1:3) = OFCmat(1:3,1:3)
        else if (FI == 'I') then
          sg%T(1:3,1:3) = OBCmat(1:3,1:3)
        else if (FI == 'A') then
          sg%T(1:3,1:3) = OBaseAmat(1:3,1:3)
        else
          write(*,*) 'FI  is '//FI
          write(*,'(A)') 'err: It must be C, F, I, or A for orthorhombic'
          stop 1
        endif
      else if(SGn >= 75 .and. SGn <= 142) then
        if(trim(ext) /= 'stdbox') then
          write(*,'(A)') 'err: not stdbox.'
          stop 1
        endif
        if(FI == 'I') then
          sg%T(1:3,1:3) = BCTmat(1:3,1:3)
        else
          write(*,*) 'Bad SGn = ',SGn
          write(*,'(A)') 'err: crash.. FC should be I only.'
          stop 1
        endif
      else if(SGn >= 143 .and. SGn <= 167) then
        if(trim(ext) == 'stdbox') then
          write(*,'(A)') 'err: trigonal cannot be stdbox. It must be mat6 (Y) or mat2 (inverted Y)'
          stop 1
        endif
        write(*,*) 'triginal: SGn= ',SGn
        if(FI == 'R' .and. trim(ext) == 'mat2') then
          sg%T(1:3,1:3) = Mat2Trigonal(1:3,1:3)
        else if(FI == 'R' .and. trim(ext) == 'mat6') then
          sg%T(1:3,1:3) = Mat6Trigonal(1:3,1:3)
        else
          write(*,'(A)') 'err: crash, should be R only.'
          stop 1
        endif
      else if(SGn >= 168 .and. SGn <= 194) then
        write(*,*) 'Hexagonal case does not a chance to use stdbox'
        write(*,'(A)') 'Hexagonal'
      else if(SGn >= 195 .and. SGn <= 230) then
        write(*,*) 'cubic: SGn= ',SGn
        if(FI == 'F') then
          sg%T(1:3,1:3) = FCCmat(1:3,1:3)
        else if(FI == 'I') then
          sg%T(1:3,1:3) = BCCmat(1:3,1:3)
        else
          write(*,'(A)') 'err: crash.. should be F or I only.'
          stop 1
        endif
      else
         write(*,*) 'Bad SGn = ',SGn
         write(*,'(A)') 'err: SGn problem'
         stop 1
      endif
      sg%origvec(1:3) = zero
    else
      open(unit=dtu,file=trim(qtablepath)//'/'//trim(dtable),status='old',action='read')
      read(dtu,*)
      found=.false.
      do
        if(found) exit
        read(dtu,DSF) tmp1

        if(trim(tmp1) /= '') then
          s = scan(tmp1,' ')
          sgkeyword = tmp1(1:s-1)
          tmp2 = tmp1(s+1:DSL)
          s = scan(tmp2,' ')
          sglabel = tmp2(1:s-1)
        endif

        if(trim(sgkeyword) == 'END') then
          exit
        endif
        if(trim(sglabel) == trim(targetlabel)) then
          found = .true.
        else

        endif
      enddo

      if(.not. found) then
        write(*,'(A)') 'SG under consideration is '//trim(targetlabel)
        write(*,'(A)') 'err: in update_derivsg_table, cannot find SG label.'
        stop 1
      endif

      read(dtu,*) sg%T(1,1:3)
      read(dtu,*) sg%T(2,1:3)
      read(dtu,*) sg%T(3,1:3)
      write(*,*) 'To find keyword origvec:'
      read(dtu,DSF) line

      read(line,*) tmp1
      if(trim(tmp1) == '#optional') then
        read(line,*) tmp1,tmp2
        write(*,'(A)') 'tmp1,tmp2='//trim(tmp1),trim(tmp2)
        if(trim(tmp2) == 'origvec') then
          write(*,*) 'keyword origvec found.'
          read(dtu,*) sg%origvec(1:3)
        else
          write(*,*) 'I: keyword origvec not found'
          sg%origvec(1:3) = zero
        endif
      else
        write(*,*) 'II: keyword origvec not found'
        sg%origvec(1:3) = zero
      endif
      write(*,'(A,3F20.10)') 'sg%origvec(1:3)=', sg%origvec(1:3)
    endif

    sg%invT = inv3x3(sg%T(1,1))

    invT(1:3,1:3) = sg%invT(1:3,1:3)

    write(*,'(A)') 'Prepare to transform the translational and rotation parts: T='
    write(*,'(3F12.5)') sg%T(1,1:3)
    write(*,'(3F12.5)') sg%T(2,1:3)
    write(*,'(3F12.5)') sg%T(3,1:3)
    write(*,'(A,F12.5)') 'det3(sg%T) = ',det3(sg%T)
    write(*,'(A)') 'Prepare to transform the translational and rotation parts: invT=T^{-1}'
    write(*,'(3F12.5)') invT(1,1:3)
    write(*,'(3F12.5)') invT(2,1:3)
    write(*,'(3F12.5)') invT(3,1:3)

    nop = sg%nop
    write(*,*) 'nop = ',nop
    do i = 1, nop

      blk4x4(1:4,1:4) = sg%op(1:4,1:4,i)
      t1(1:3) = blk4x4(2:4,1)
      R1(1:3,1:3) = blk4x4(2:4,2:4)

      R2 = BAinvB(invT(1,1),R1(1,1))

      origvec(1:3) = sg%origvec(1:3)
      call matvecn(3,R1(1,1),origvec(1),t3(1))

      t1(1:3) = t1(1:3) + t3(1:3) - origvec(1:3)

      call matvecn(3,invT(1,1),t1(1),t2(1))

      call foldtofirstzone3(t2(1))

      newblk4x4(1:4,1:4) = zero
      newblk4x4(1,1) = one
      newblk4x4(2:4,1) = t2(1:3)
      newblk4x4(2:4,2:4) = R2(1:3,1:3)

      write(*,'(A,I5)') 'Pre || Post:  i =',i
      write(*,'(4F12.5,A,4F12.5)') blk4x4(2,1:4),' || ', newblk4x4(2,1:4)
      write(*,'(4F12.5,A,4F12.5)') blk4x4(3,1:4),' || ', newblk4x4(3,1:4)
      write(*,'(4F12.5,A,4F12.5)') blk4x4(4,1:4),' || ', newblk4x4(4,1:4)

      sg%op(1:4,1:4,i) = newblk4x4(1:4,1:4)
    enddo

    write(*,*)
    filename = 'echo-ext-operations.dat'
    write(*,'(A)') 'To produce '//trim(filename)

    pru = newunit()
    open(unit=pru,file=trim(filename),status='replace')

    write(pru,'(A)') '---------------------------------------'
    write(pru,'(A)') '    SG = '//trim(targetlabel)
    write(*,'(A)') '    SG = '//trim(targetlabel)
    write(pru,'(A)') '    Y = y'
    write(pru,'(A)') '    SGbase(Y,SG)%sglabel    ='''//trim(targetlabel)//''''
    write(pru,'(A)') '    SGbase(Y,SG)%PG         ='''//trim(sg%pg)//''''
    write(pru,'(A)') '    SGbase(Y,SG)%ncentering ='//trim(N2str(sg%ncentering))
    do i = 1, sg%ncentering
      do j = 1, 3
        f = sg%centering(j,i)

        call frac2nd(f,nf,df)
        if(abs(df-one) < 1.0d-6) then

          string(j)=trim(N2str(nint(nf)))//'.0d0'
        else
          string(j)=trim(N2str(nint(nf)))//'.0d0/'//trim(N2str(nint(df)))//'.0d0'
        endif
      enddo
      write(pru,'(A,A16,A,A16,A,A16,A)') '    SGbase(Y,SG)%centering(1:3,'//trim(N2str(i))//') = (/',trim(string(1)),',',trim(string(2)),',',trim(string(3)),' /)'
    enddo
    write(pru,'(A)') '    SGbase(Y,SG)%PGNumOper  ='//trim(N2str(sg%PGNumOper))
    write(pru,'(A)') '    SGbase(Y,SG)%symmorphic ='//trim(N2str(sg%symmorphic))
    if(sg%symmorphic == 0) then
      do i = 1, sg%PGNumOper
        do j = 1, 3
          f = sg%op(j+1,1,i)

          call frac2nd(f,nf,df)
          if(abs(df-one) < 1.0d-6) then

            string(j)=trim(N2str(nint(nf)))//'.0d0'
          else
            string(j)=trim(N2str(nint(nf)))//'.0d0/'//trim(N2str(nint(df)))//'.0d0'
          endif
        enddo
        write(pru,'(A,A16,A,A16,A,A16,A)') '    SGbase(Y,SG)%nsy(1:3,'//trim(N2str(i))//') = (/',trim(string(1)),',',trim(string(2)),',',trim(string(3)),' /)'
      enddo
    endif
    write(pru,'(A)') '    SGbase(Y,SG)%nop        ='//trim(N2str(sg%nop))
    write(pru,'(A)') '    SGbase(Y,SG)%assigned   =.TRUE.'
    write(pru,*)
    close(pru)
    write(*,'(A)') trim(filename)//' is generated.'
    write(*,*)

    write(*,*) 'original sg%sglabel = '//trim(sg%sglabel)
    sg%sglabel = trim(targetlabel)
    write(*,'(A)') 'Updated: sg%sglabel = '//trim(sg%sglabel)
  end subroutine update_derivsg_table

  subroutine assign_m01(m01,ns0,chPOSCAR,pPOSCAR,distance_thres)
    integer :: ns0,m01(ns0)
    type(supercell) :: chPOSCAR,pPOSCAR
    integer :: i,j,ii,jj,kk
    real(double) :: v1(3),v2(3),v3(3),newf(3),distance_thres
    m01(:) = -1
    do i = 1, ns0
      call frac2abs(chposcar%a,chposcar%at(i)%f,v1)
      do j = 1, pPOSCAR%n
        do kk = -1, 1
          do jj = -1, 1
            do ii = -1, 1
              newf = pPOSCAR%at(j)%f + (/ii,jj,kk/)
              call frac2abs(pPOSCAR%a,newf,v2)
              v3 = v2-v1
              if(vecmag3(v3) < distance_thres) then
                m01(i) = j
                exit
              endif
            enddo
          enddo
        enddo
      enddo
      if(m01(i) == -1) then
        write(*,*) 'atom', i, ' (in ch.poscar) cannot be found in pPOSCAR.'
        write(*,'(A)') 'err: serious mapping error in assign_m01'
        stop 1
      endif
    enddo
  end subroutine assign_m01

  subroutine assign_m12_tvec(tvec,m12,ns1,s1,s2,distance_thres)
    real(double) :: tvec(3)
    integer :: ns1,m12(ns1)
    type(supercell) :: s1,s2
    real(double) ::distance_thres
    integer :: i,j,ii,jj,kk
    real(double) :: v1(3),v2(3),v3(3),newf(3)
    real(double) :: mindis,dis
    m12(1:ns1) = -1

    do i = 1, ns1
      call frac2abs(s1%a(1,1),s1%at(i)%f(1),v1(1))

      v1(1:3) = v1(1:3) + tvec(1:3)

      mindis = 1d100
      do j = 1, s2%n

        do kk = -1, 1
          do jj = -1, 1
            do ii = -1, 1
              newf = s2%at(j)%f + (/ii,jj,kk/)
              call frac2abs(s2%a(1,1),newf(1),v2(1))
              v3 = v2-v1
              dis = vecmag3(v3(1))
              if(dis < mindis) then
                mindis = dis
              endif
              if(dis < distance_thres) then
                m12(i) = j
                exit
              endif
            enddo
          enddo
        enddo
      enddo

      if(m12(i) == -1) then
        write(*,*) 'mindis = ',mindis
        write(*,*) 'atm index=', i, ', cannot be found in supercell.poscar.'
        write(*,'(A)') 'err: serious mapping error in assign_m12'
        stop 1
      endif
    enddo

  end subroutine assign_m12_tvec

  subroutine assign_m21(s2,s1,m21,ns2)
    type(supercell) :: s1,s2
    integer :: ns2,ns1
    integer :: i,j,m21(ns2)
    real(double) :: v2(3),v1(3),diffv(3),a1(3,3),inva1(3,3)
    real(double) :: sol(3),devia(3)

    m21(1:ns2) = -1

    a1(1:3,1:3) = s1%a(1:3,1:3)

    inva1 = inv3x3(a1)

    ns1 = s1%n

    do i = 1, ns2
      call frac2abs(s2%a,s2%at(i)%f,v2)
      do j = 1, ns1
        call frac2abs(s1%a,s1%at(j)%f,v1)

        diffv = v2-v1
        call matvecn(3,inva1,diffv,sol)

        devia(1:3) = abs(sol(1:3)-nint(sol(1:3)))

        if(vecmag3(devia(1)) < 1.d-10) then
          m21(i) = j

        endif
      enddo
      if(m21(i) == -1) then
        write(*,*) 'i = ',i
        write(*,*) 'cannot find an equivalent atom in the pPOSCAR'
        write(*,'(A)') 'err: serious translational error.'
        stop 1
      endif
    enddo
  end subroutine assign_m21

  subroutine get_newrot(newrot,pposcar,sp,op)
    type(onesg) :: sp
    integer :: op
    real(double) :: M1(3,3),newrot(3,3),rot(3,3),lattmat(3,3),invlattmat(3,3)
    type(supercell) :: pposcar
    Rot(1:3,1) = sp%op(2:4,2,Op)
    Rot(1:3,2) = sp%op(2:4,3,Op)
    Rot(1:3,3) = sp%op(2:4,4,Op)

    LattMat(1:3,1:3) = pPOSCAR%a(1:3,1:3)

    InvLattMat = inv3x3(LattMat(1,1))
    call matmatn(3,LattMat,Rot,M1)
    call matmatn(3,M1,InvLattMat,NewRot)
  end subroutine get_newrot

  subroutine get_fundmental_const_in_internal_units(MLTIQ,bu,fc)
    type(internalunitinsi) :: bu
    type(internalfundamentalconstant) :: fc
    real(double) :: conv
    real(double) :: MLTIQ(5)

    bu%M = MLTIQ(1)
    bu%L = MLTIQ(2)
    bu%T = MLTIQ(3)
    bu%I = MLTIQ(4)
    bu%Q = MLTIQ(5)

    conv = (bu%M)*(bu%L)/(bu%I**2)/(bu%T**2)
    fc%internalmu0 = mu0/conv

    conv = bu%I*bu%T
    fc%internalecharge = EChg/conv

    conv = bu%I*bu%L**2
    fc%internalmuB = BohrMagneton/conv

    conv = bu%M*bu%L**2/bu%T
    fc%internalhbar = hbar/conv

    conv = bu%M*(bu%L)**2d0/(bu%T**2d0)/bu%Q
    fc%internalkB = kBoltz/conv
  end subroutine get_fundmental_const_in_internal_units

  subroutine cyclic_index(id,N,nproc,b,cycind,ns)
    integer :: id,N,nproc,b
    integer :: cycind(N)
    integer :: nb,k,copyn,ind,j,ns

    if(id < 1 .or. id > nproc) then
      write(*,*) 'invalid processor number : id =',id
      write(*,'(A)') 'err: check you processor ID'
      stop 1
    endif
    if(nproc < 1) then
      write(*,*) 'invalid number of processors : nproc = ',nproc
      write(*,'(A)') 'err: check the number of processors'
      stop 1
    endif
    if(b < 1) then
      write(*,*) 'invalid block number : b =',b
      write(*,'(A)') 'err: check the block size'
      stop 1
    endif

    nb = N/b
    if(nb*b /= N) then
      nb = nb + 1
    endif

    ns = 0
    do j = 1, nb
      ind = mod(j-1,nproc)+1
      if(id == ind) then
        if(j == nb) then

          copyn = N - (nb-1)*b
          do k = 1, copyn
            ns = ns + 1
            cycind(ns) = (j-1)*b + k
          enddo
        else
          do k = 1, b
            ns = ns + 1
            cycind(ns) = (j-1)*b + k
          enddo
        endif
      endif
    enddo
  end subroutine cyclic_index

  subroutine get_grid_dim(nproc,ngr,ngc)
    integer :: nproc,ngr,ngc
    integer :: nstart
    logical :: done
    integer :: nrow,nscan

    if(nproc < 1) then
       write(*,'(A)') 'err: nproc must be greater than 0'
       stop 1
    endif
    nstart = nint(sqrt(nproc*one))
    done = .false.
    nscan = nstart
    do
      if(done) exit
      nrow = nproc/nscan
      if(nscan*nrow == nproc) then
        done = .true.
        ngr = nrow
        ngc = nscan
      else
        nscan = nscan - 1
      endif
    enddo
  end subroutine get_grid_dim

  function D1R_assign(alpha,beta,gamm) result (r)
    complex(double) :: r(3,3)
    real(double) :: alpha,beta,gamm
    r(1,1) = (1d0/2d0)*(1d0+cos(beta))*exp(-Imag*(1*alpha + 1*gamm))
    r(1,2) = (-1d0/sqrt2)*(sin(beta))*exp(-Imag*(1*alpha + 0*gamm))
    r(1,3) = (1d0/2d0)*(1d0-cos(beta))*exp(-Imag*(1*alpha + (-1)*gamm))
    r(2,1) = (1d0/sqrt2)*(sin(beta))*exp(-Imag*(0*alpha + 1*gamm))
    r(2,2) = (1d0/1d0)*(cos(beta))*exp(-Imag*(0*alpha + 0*gamm))
    r(2,3) = (-1d0/sqrt2)*(sin(beta))*exp(-Imag*(0*alpha + (-1)*gamm))
    r(3,1) = (1d0/2d0)*(1d0-cos(beta))*exp(-Imag*((-1)*alpha + (1)*gamm))
    r(3,2) = (1d0/sqrt2)*(sin(beta))*exp(-Imag*((-1)*alpha + (0)*gamm))
    r(3,3) = (1d0/2d0)*(1d0+cos(beta))*exp(-Imag*((-1)*alpha + (-1)*gamm))
  end function D1R_assign

  function yl1m(theta,phi) result (s)
    complex(double) :: s(3)
    real(double) :: theta,phi

    s(1) = -(3d0/(8d0*pi))**0.5d0*sin(theta)*exp(Imag*phi)
    s(2) = (3d0/(4d0*pi))**0.5d0 * cos(theta)
    s(3) = (3d0/(8d0*pi))**0.5d0*sin(theta)*exp(-Imag*phi)
  end function yl1m

  function fermifunc(E,mu,T)
    real(double) :: E,mu,T,fermifunc,u
    u = (E-mu)/(kBoltzDeV*T)
    if(u < zero) then
      fermifunc = 1d0/(exp(u) + 1d0)
    else
      fermifunc = exp(-u)/(1d0 + exp(-u))
    endif
  end function fermifunc

  function dfde(E,mu,T)
    real(double) :: dfde,E,mu,T,u

    u = (E-mu)/(kBoltzDeV*T)

    if(u < zero) then
      dfde = -(1d0/(kBoltzDeV*T))*exp(u)/(1d0+exp(u))**2d0
    else
      dfde = -(1d0/(kBoltzDeV*T))*exp(-u)/(1d0+exp(-u))**2d0
    endif
  end function dfde

  function entropy_func(omega,T)
    real(double) :: entropy_func,omega,s,u,T

    if(abs(omega/(1.0d12*2*pi)) < 1.0d-2) then
      entropy_func = 0.0d0
      return
    endif

    s = hbar*omega/(kBoltz*T)
    u = exp(s)
    entropy_func = s/(u-1) - log(1-(one/u))
  end function entropy_func

  function fenergy_func(omega,T)
    real(double) :: fenergy_func,omega,T

    if(abs(omega/(1.0d12*2*pi)) < 2.0d-2) then

      fenergy_func = 0.0d0
      return
    endif

    fenergy_func = kBoltz*T*log(2*sinh(hbar*omega/(2*kBoltz*T)))
    fenergy_func = fenergy_func*J2eV
  end function fenergy_func

  function Planck_nfunc(freq,T)
    real(double) :: freq,T,Planck_nfunc,r
    r = hPlanck*freq*invcm/(kBoltz*T)

    Planck_nfunc = 0.5d0/tanh(0.5d0*r) - 0.5d0
  end function Planck_nfunc

  subroutine fillall_rhombopara(rstruc)
    type(rhombopara) :: rstruc
    real(double) :: ah,ch,r,sinalpha,cosalpha,alpha,alphadeg,a,s,d,h,t
    real(double) :: rot(3,3),rot2(3,3)

    ah = rstruc%ah
    ch = rstruc%ch
    r = ch/ah

    rstruc%caratio = r

    cosalpha = (2d0*r**2d0-3d0)/(2d0*r**2d0+6d0)
    rstruc%cosalpha = cosalpha

    sinalpha = 3d0*sqrt(4d0*r**2d0 + 3d0)/(2d0*r**2d0 + 6d0)

    alpha = atan2(sinalpha,cosalpha)
    rstruc%alpha = alpha

    alphadeg = alpha*rad2deg
    rstruc%alphardeg = alphadeg

    a = (ah/3d0)*sqrt( r**2d0+3d0   )
    rstruc%ar = a

    s = cosalpha
    d = a*sin(alpha/2d0)
    h = 2d0*d/sqrt3
    t = a*sqrt( (1d0 + 2d0*s)/3d0)

    rstruc%mat0(1:3,1) = (/ah,zero,zero/)
    rstruc%mat0(1:3,2) = (/-ah/two,ah*sqrt3/two,zero/)
    rstruc%mat0(1:3,3) = (/zero,zero,ch/)

    rstruc%mat1(1:3,1) = (/zero,h,t/)
    rstruc%mat1(1:3,2) = (/-d,-h/2d0,t/)
    rstruc%mat1(1:3,3) = (/d,-h/2d0,t/)

    rot = neta_2_rotm( (/zero,zero,one/), -120.0d0)

    call matmatn(3, rot(1,1),rstruc%mat1(1,1),rstruc%mat2(1,1))

    call matmatn(3, rot(1,1),rstruc%mat2(1,1),rstruc%mat3(1,1))

    rot2 = neta_2_rotm( (/zero,zero,one/), 60.0d0)
    call matmatn(3, rot2(1,1),rstruc%mat3(1,1),rstruc%mat4(1,1))

    call matmatn(3, rot(1,1),rstruc%mat4(1,1),rstruc%mat5(1,1))

    call matmatn(3, rot(1,1),rstruc%mat5(1,1),rstruc%mat6(1,1))

  end subroutine fillall_rhombopara

  subroutine print_all_rhombopara(r)
    type(rhombopara) :: r

    write(*,*) 'A comprehensive list of rhombohedral parameters:'
    write(*,*) 'ah = ',r%ah
    write(*,*) 'ch = ',r%ch
    write(*,*) 'eta = ch/ah = ',r%caratio
    write(*,*) 'ar = ',r%ar
    write(*,*) 'alpha (deg) = ',r%alphardeg
    write(*,*) 'alphar (radian) = ',r%alpha
    write(*,*) 'cosalpha = ',r%cosalpha
    write(*,*) 'd = ',r%ah/two
    write(*,*) 'h = ',r%ah/sqrt3
    write(*,*) 't = ',r%ch/three
    write(*,*) 'eta for BS = ',(one + 4*cos(r%alpha))/(2+ 4*cos(r%alpha))
    write(*,*) 'volume = ',  r%ah**2d0*r%ch/sqrt(12.0d0)
    write(*,'(A)') '  matrix0: hexaganal'
    write(*,'(3F18.12)') r%mat0(1:3,1)
    write(*,'(3F18.12)') r%mat0(1:3,2)
    write(*,'(3F18.12)') r%mat0(1:3,3)

    write(*,'(A)') '  matrix2: inverted Y (PWSCF)'
    write(*,'(3F18.12)') r%mat2(1:3,1)
    write(*,'(3F18.12)') r%mat2(1:3,2)
    write(*,'(3F18.12)') r%mat2(1:3,3)

    write(*,'(A)') '  matrix6: upright Y (default for trigonal/rhombohedral)'
    write(*,'(3F18.12)') r%mat6(1:3,1)
    write(*,'(3F18.12)') r%mat6(1:3,2)
    write(*,'(3F18.12)') r%mat6(1:3,3)
  end subroutine print_all_rhombopara

  subroutine internal_ah_ch_rhombopara(ah,ch,rstruc)
    real(double) :: ah,ch
    type(rhombopara) :: rstruc
    rstruc%ah = ah
    rstruc%ch = ch
    call fillall_rhombopara(rstruc)
    call print_all_rhombopara(rstruc)
  end subroutine internal_ah_ch_rhombopara

  subroutine ah_ch_rhombopara(ah,ch,rstruc)
    real(double) :: ah,ch
    type(rhombopara) :: rstruc
    call internal_ah_ch_rhombopara(ah,ch,rstruc)
  end subroutine ah_ch_rhombopara

  subroutine a_alphadeg_rhombopara(a,alphadeg,rstruc)
    real(double) :: a,alphadeg
    type(rhombopara) :: rstruc
    real(double) :: ah,ch

    ah = 2d0*a*sin(alphadeg*deg2rad/2.0d0)
    ch = a*sqrt(3.0d0*(1.0d0+2.0d0*cos(alphadeg*deg2rad)))
    call internal_ah_ch_rhombopara(ah,ch,rstruc)
  end subroutine a_alphadeg_rhombopara

  subroutine d_and_t_rhombopara(d,t,rstruc)
    type(rhombopara) :: rstruc
    real(double) :: d,t
    real(double) :: ah,ch

    ah = two*d
    ch = three*t
    call internal_ah_ch_rhombopara(ah,ch,rstruc)
  end subroutine d_and_t_rhombopara

  subroutine least_sq_fit(ndeg,ndata,offset,xval,yval,s,M,MT,MTM,MTy,sol,ipiv)
    integer :: ndeg,ndata,ndeg1
    real(double) :: offset
    real(double) :: xval(ndata),yval(ndata)
    real(double) :: M(ndata,0:ndeg),MT(0:ndeg,ndata),MTM(0:ndeg,0:ndeg),MTy(0:ndeg)
    real(double) :: sol(0:ndeg,1)
    real(double) :: s(0:ndeg)
    integer :: ipiv(0:ndeg)
    integer :: i,j,info
    real(double) :: alpha,beta

    do i = 1, ndata
      do j = 0, ndeg
        if(j == 0) then
          M(i,j) = 1.0D0
        else
          M(i,j) = (xval(i)-offset)**(j*1.0D0)
        endif
      enddo
    enddo

    do i = 0, ndeg
      do j = 1, ndata
        MT(i,j) = M(j,i)
      enddo
    enddo
    ndeg1=ndeg+1
    alpha = One
    beta = Zero
    call dgemm('N','N',ndeg1,ndeg1,ndata,1.0D0,MT(0,1),ndeg1,M,ndata,0.0D0,MTM(0,0),ndeg1)
    alpha = One
    beta = Zero
    call dgemv('N',ndeg1,ndata,alpha,MT(0,1),ndeg1,yval(1),1,beta,MTy(0),1)
    sol(:,1) = MTy(:)

    CALL dgesv(ndeg1,1,MTM(0,0),ndeg1,ipiv(0),sol(0,1),ndeg1,info)
    if(info /= 0) then
      write(*,*) 'In least_sq_fit, info = ',info
      write(*,'(A)') 'err: info not zero.'
      stop 1
    endif
    s(:) = sol(:,1)
    do i = 0, ndeg

    enddo
  end subroutine least_sq_fit

  subroutine fx_fpx_least_sq_fit(x,a,ndeg,offset,fx,fpx)
    real(double) :: fx,fpx,x,a(0:ndeg),offset
    integer :: i,ndeg
    fx = zero
    fpx = zero
    do i = 0, ndeg
      fx = fx + a(i)*(x-offset)**(i*one)
      if(i > 0) then
        fpx = fpx + i*a(i)*(x-offset)**(i-one)
      endif
    enddo
  end subroutine fx_fpx_least_sq_fit

  function strainvolratio(e)
    real(double) :: e(6),strainvolratio
    real(double) :: M(3,3)
    M(1:3,1) = (/ one + e(1), e(6)/two, e(5)/two /)
    M(1:3,2) = (/ e(6)/two, one + e(2), e(4)/two /)
    M(1:3,3) = (/ e(5)/two, e(4)/two, one + e(3) /)
    strainvolratio = det3(M(1,1))
  end function strainvolratio

  subroutine create_primitive_cell(sp,s,t,kstar_arr)
    type(onesg) :: sp
    type(supercell) :: s
    type(supercell) :: t
    integer :: i,j,k
    real(double) :: vec(4,1),prod(4,1),p(3),absp(3)
    logical :: newatom
    integer :: indexu
    type(onekstar) :: kstar_arr(s%n)
    integer :: runindex,opind
    real(double) :: v1(4),v2(4),v3(3),matop(4,4),p1(3),p2(3),newf(3),rotm(3,3)

    t%n = 0
    t%a = s%a
    call real2recip(t%a,t%b)

    indexu = newunit()

    open(unit=indexu,file='ineq_atms_loc.dat',status='replace')

    t%nsp = s%nsp
    t%zarr(1:t%nsp) = s%zarr(1:s%nsp)

    do i = 1, s%n
      k = 0
      do j = 1, sp%nop

        vec(1,1) = 1d0
        vec(2:4,1) = s%at(i)%f

        call matmatmnp(4,4,1,sp%op(1:4,1:4,j),vec(1:4,1),prod(1:4,1))

        p = prod(2:4,1)
        call frac2abs(s%a,p,absp)
        call tryput(s%at(i)%z,absp,t,newatom)
        if(newatom) then
          k = k + 1
          if(k == 1) then

            write(indexu,'(I8)') t%n
          else
            write(indexu,'(A8)') 'equiv'
          endif
          kstar_arr(i)%op(k) = j
          kstar_arr(i)%ind1(k) = t%n
          kstar_arr(i)%trans(1:3,k) = (/0,0,0/)
        endif
      enddo
      kstar_arr(i)%n = k
    enddo
    close(unit=indexu)
    write(*,*)

    write(*,'(A,I5)') 'Full poscar has totaln = ',t%n

    write(*,*)

    runindex = 0
    do i = 1, s%n
      v1(1:4) = (/one,s%at(i)%f(1:3)/)
      do j = 1, kstar_arr(i)%n
        runindex = runindex + 1
        opind = kstar_arr(i)%op(j)
        matop(1:4,1:4) = sp%op(1:4,1:4,opind)
        call matvecn(4,matop(1,1),v1(1),v2(1))

        v3(1:3) = v2(2:4)
        call foldtofirstzone3(v3(1))

        call frac2abs(s%a,v3(1),p1(1))

        call foldtofirstzone3(t%at(runindex)%f)

        call frac2abs(t%a,t%at(runindex)%f,p2(1))

        if(vecmag3(p1-p2) > 1.0d-5) then
          write(*,'(A,I5,6F20.10,F20.10)') 'runindex,p1,p2=',runindex,p1(1:3),p2(1:3),vecmag3(p1-p2)
          write(*,'(A)') 'err: star of k is wrong?'
          stop 1
        endif

        rotm(1:3,1:3) = matop(2:4,2:4)
        call matvecn(3,rotm(1,1),s%at(i)%force(1),newf(1))
        t%at(runindex)%force(1:3) = newf(1:3)
      enddo
    enddo
  end subroutine create_primitive_cell

  subroutine find_group_of_k(ns0,chposcar,pposcar,sg,m01,kgroup_arr)
    integer :: ns0
    type(supercell) :: chposcar,pposcar
    type(onesg) :: sg
    integer :: m01(ns0)
    type(onekgroup) :: kgroup_arr(ns0)
    integer :: i,j,gs,ii,jj,kk
    integer :: ind1
    real(double) :: fracA(3),vec(4,1),prod(4,1),p(3),q(3),diffv(3),dis
    logical :: found

    integer :: absnrange=5

    ns0 = chposcar%n
    do i = 1, ns0
      ind1 = m01(i)
      gs = 0
      fracA = pposcar%at(ind1)%f
      do j = 1, sg%nop
        vec(1,1) = one
        vec(2:4,1) = fracA
        call matmatmnp(4,4,1,sg%op(1:4,1:4,j),vec(1:4,1),prod(1:4,1))
        p = prod(2:4,1)
        found = .false.

        do ii = -absnrange, absnrange
          do jj = -absnrange, absnrange
            do kk = -absnrange, absnrange

              q = p + (/ii,jj,kk/)
              diffv = q-fracA
              dis = vecmag3(diffv)
              if(dis < 1d-8) then
                found = .true.
                gs = gs + 1
                kgroup_arr(i)%op(gs) = j
              endif
            enddo
          enddo
        enddo
      enddo
      kgroup_arr(i)%n = gs
    enddo
  end subroutine find_group_of_k

  function check_inv(sg)
    integer :: check_inv
    type(onesg) :: sg
    integer :: i,storedi
    real(double) :: invm(3,3),diffm(3,3),tau(3),m(3,3)

    check_inv = 0
    invm(1:3,1:3) = zero
    invm(1,1) = -one
    invm(2,2) = -one
    invm(3,3) = -one
    do i = 1, sg%nop
      m(1:3,1:3) = sg%op(2:4,2:4,i)
      diffm = m(1:3,1:3) - invm(1:3,1:3)
      tau(1:3) = sg%op(2:4,1,i)
      if(matmagn(3,diffm(1,1)) < 1.0d-8 .and. vecmag3(tau(1)) < 1.0d-8) then

        storedi = i
        check_inv = 1
      endif
    enddo
    write(*,*)
    if(check_inv == 0) then
      write(*,'(A)') 'The space group does not have an inversion symmetry'
    else
      write(*,'(A)') 'The space group has an inversion symmetry:'
      write(*,'(A,I4)') 'Operation number is ',storedi
      write(*,'(A)') 'Operation matrix is:'
      write(*,'(A,9F20.10)') 'm(1,1:3) = ',sg%op(2,2:4,storedi)
      write(*,'(A,9F20.10)') 'm(2,1:3) = ',sg%op(3,2:4,storedi)
      write(*,'(A,9F20.10)') 'm(3,1:3) = ',sg%op(4,2:4,storedi)
      write(*,'(A,3F20.10)') 'tau(1:3) = ',sg%op(2:4,1,storedi)
    endif
  end function check_inv

  subroutine get_supercell_atom_map(u,pposcar,superc,map)

    real(double) :: u(4,4)
    type(supercell) :: pposcar,superc,tmp1
    integer :: map(superc%n)
    integer :: j,k
    real(double) :: jabsc(3),jfracc(3), prod(4,4),v(4,1)
    real(double) :: dev,testt(3),diffv(3),rundev
    integer :: n1scan,n2scan,n3scan
    real(double) :: summap,summap_analytic
    integer :: ns2
    integer,allocatable :: invmap(:)

    tmp1%a = superc%a
    call real2recip(tmp1%a,tmp1%b)
    tmp1%n = superc%n
    do j = 1, tmp1%n
      tmp1%at(j)%z = superc%at(j)%z
      tmp1%at(j)%f = superc%at(j)%f
    enddo

    do j = 1, tmp1%n

      call frac2abs(tmp1%a, tmp1%at(j)%f(1:3),jabsc)
      call abs2frac(pposcar%b,jabsc,jfracc)

      v(1,1) = one
      v(2:4,1) = jfracc(1:3)

      call matmatmnp(4,4,1,u(1:4,1:4),v(1:4,1),prod(1:4,1))

      jfracc(1:3) = prod(2:4,1)
      call frac2abs(pposcar%a,jfracc,jabsc)
      call abs2frac(tmp1%b,jabsc,tmp1%at(j)%f)

      call foldtofirstzone3(tmp1%at(j)%f)
    enddo

    allocate(invmap(superc%n))

    map(:) = 0
    invmap(:) = 0

    do j = 1, tmp1%n

      do k = 1, superc%n

        dev = 1d100
        do n1scan = -1, 1
          do n2scan = -1, 1
            do n3scan = -1, 1
              testt = tmp1%at(j)%f + (/n1scan,n2scan,n3scan/)
              diffv = testt - superc%at(k)%f
              rundev = vecmag3(diffv)
              if(rundev < dev) then
                dev = rundev
              endif
            enddo
          enddo
        enddo

        if(dev < dev_tol_get_supercell_atom_map) then
          if(map(j) /= 0) then
            write(*,*) 'Two atoms in the old list match with this atom j'
            write(*,'(A)') 'err: Check again.'
            stop 1
          endif
          map(j) = k
          if(invmap(k) /= 0) then
            write(*,*) 'j,k = ',j,k
            write(*,*) 'k is pre-mapped to a previous atom=',invmap(k)
            write(*,'(A)') 'err: Check again.'
            stop 1
          endif
          invmap(k) = j
        endif
      enddo
      if(map(j) == 0) then
        write(*,*) 'problematic j = ',j
        write(*,*) 'something is wrong. mapj must be nonzero'
        write(*,'(A)') 'err: stop here.'
        stop 1
      else

      endif
    enddo

    deallocate(invmap)

    ns2 = superc%n
    summap = 0.0d0
    do j = 1, ns2
      if(map(j) < 0 .or. map(j) > ns2) then
        write(*,*) 'j,map(j)=',j,map(j)
        write(*,'(A)') 'err: sanity check fails.'
        stop 1
      endif
      summap = summap + map(j)
    enddo
    summap_analytic=ns2*(one+ns2)/two
    if(abs(summap - summap_analytic) > 1.0d-10) then
      write(*,*) 'summap,summap_analytic=',summap,summap_analytic
      write(*,*) 'map array sanity check fails.'
      write(*,'(A)') 'err: check.'
      stop 1
    endif
  end subroutine get_supercell_atom_map

  subroutine get_normalvec(P,normal)
    real(double) :: P(3,3),normal(4),m(3,3)
    real(double) :: a,b,c,d,sq

    m(1,1:3) = (/  P(2,1), P(3,1), one /)
    m(2,1:3) = (/  P(2,2), P(3,2), one /)
    m(3,1:3) = (/  P(2,3), P(3,3), one /)
    a = det3(m(1,1))

    m(1,1:3) = (/  P(3,1), P(1,1), one /)
    m(2,1:3) = (/  P(3,2), P(1,2), one /)
    m(3,1:3) = (/  P(3,3), P(1,3), one /)
    b = det3(m(1,1))

    m(1,1:3) = (/  P(1,1), P(2,1), one /)
    m(2,1:3) = (/  P(1,2), P(2,2), one /)
    m(3,1:3) = (/  P(1,3), P(2,3), one /)
    c = det3(m(1,1))

    m(1,1:3) = (/  P(1,1), P(2,1), P(3,1) /)
    m(2,1:3) = (/  P(1,2), P(2,2), P(3,2) /)
    m(3,1:3) = (/  P(1,3), P(2,3), P(3,3) /)
    d = -det3(m(1,1))
    write(*,'(A,4F15.8)') 'a,b,c,d=',a,b,c,d
    sq = sqrt(a**2 + b**2 + c**2)
    if(sq < 1.d-10) then
      write(*,*) 'sq =',sq
      write(*,*) 'A,B,C are on the same line?'
      write(*,'(A)') 'err: Check if A, B, and C are on the same line.'
      stop 1
    endif
    normal(1:4) =  (one/sq)*(/  a,b,c,d /)
  end subroutine get_normalvec

  subroutine find_ellipse_area1(zv,area,r1,r2)
    complex(double) zv(3),s(3)
    real(double) :: r1,r2
    real(double) :: magarr(3),p0,area,det
    integer :: i, j,ind(3)
    real(double) :: r(3),phase(3)
    integer,parameter :: LRWORK=9
    real(double) :: a(3,2),at(2,3),ata(2,2),atainv(2,2),q(2,3),qt(3,2),L(3,3),eig(3),RWORK(LRWORK)
    integer :: info

    do i = 1, 3
      magarr(i) = abs(zv(i))
    enddo

    do i = 1, 3
      ind(i) = i
    enddo

    call dbl_sort(3,magarr(1),ind(1),-1)
    do i = 1, 3

      s(i) = zv(ind(i))
    enddo

    p0 = atan2(  aimag(s(1)), real(s(1)))
    do i = 1, 3
      s(i) = s(i)* exp(-Imag* p0)
      r(i) = abs(s(i))
      phase(i) = atan2(  aimag(s(i)), real(s(i)))

    enddo

    a(1:3,1) = (/  r(1),  r(2) * cos(phase(2)),  r(3)*cos(phase(3)) /)
    a(1:3,2) = (/  zero, -r(2) * sin(phase(2)), -r(3)*sin(phase(3)) /)
    do i = 1, 2
      do j = 1, 3
        at(i,j) = a(j,i)
      enddo
    enddo

    call matmatmnp(2,3,2,at(1,1),a(1,1),ata(1,1))

    det = ata(1,1)*ata(2,2) - ata(1,2)*ata(2,1)
    if(abs(det) < 1.0d-14) then
      area = zero
      return
    else
      atainv(1,1) = ata(2,2)
      atainv(2,2) = ata(1,1)
      atainv(1,2) = -ata(1,2)
      atainv(2,1) = -ata(2,1)
      atainv(1:2,1:2) = atainv(1:2,1:2)/det

      call matmatmnp(2,2,3, atainv(1,1),at(1,1),q(1,1))

      do i = 1, 3
        do j = 1, 2
          qt(i,j) = q(j,i)
        enddo
      enddo

      call matmatmnp(3,2,3,qt(1,1),q(1,1),L(1,1))

      call dsyev('N','L',3,L(1,1),3,eig(1),RWORK(1),LRWORK,info)
      if(info /= 0) then
        write(*,'(A)') 'err: dsyev problem.'
        stop 1
      endif
      if(abs(eig(1)) > 1d-4) then

        write(*,*) 'Problematic zv:'
        write(*,*) 'zv(1:3)=',zv(1:3)
        write(*,*) 'a(1:3,1) = ',a(1:3,1)
        write(*,*) 'a(1:3,2) = ',a(1:3,2)
        write(*,*) 'eig(1) = ',eig(1)
        write(*,'(A)') 'err: too large magnitude. Supposed to be zero for an ellipse.'
        stop 1

      endif

      r1 = one/sqrt(eig(2))
      r2 = one/sqrt(eig(3))

      area = pi*r1*r2

    endif
  end subroutine find_ellipse_area1

  subroutine find_ellipse_area2(zv,area,r1,r2)
    complex(double) :: zv(3)
    real(double) :: area,r1,r2
    real(double) :: area2,Rz(3),Iz(3),crossp(3)

    real(double) :: a1,b1,a2,b2,a3,b3
    real(double) :: A,B,C,Ap,Cp
    real(double) :: R,S
    real(double) :: konst

    a1 = real(zv(1))
    b1 = aimag(zv(1))

    a2 = real(zv(2))
    b2 = aimag(zv(2))
    a3 = real(zv(3))
    b3 = aimag(zv(3))

    A = a1*a1 + a2*a2 + a3*a3
    B = b1*b1 + b2*b2 + b3*b3
    C = -(a1*b1 + a2*b2 + a3*b3)

    Ap = (A-B)/2d0
    Cp = C

    S = (A+B)/2d0

    R = sqrt(Ap*Ap + Cp*Cp)

    r1 = sqrt(S+R)

    if(R > S) then

      if(abs(R-S) < 1d-15) then
        r2 = zero
      else
        write(*,*) 'too large a difference.'
        write(*,'(A)') 'err: Give up.'
        stop 1
      endif
    else
      r2 = sqrt(S-R)
    endif

    konst = pi

    area = konst*r1*r2

    Rz(1:3) = real(zv(1:3))
    Iz(1:3) = aimag(zv(1:3))

    crossp = crossprod(rz(1),iz(1))
    area2 = konst*vecmag3(crossp(1))

    if( abs(area2 - area) > 1d-7) then

      write(*,*) 'area2,area=',area2,area

      write(*,'(A)') 'err: Too big difference from two appraoches of obtaining areas.'
      stop 1
    else

    endif
  end subroutine find_ellipse_area2

  subroutine read_forces(fu,abspartialfrcfile,ns1,ns2,zf,nzf)
    integer :: fu
    character(len=*) :: abspartialfrcfile
    integer :: ns1,ns2
    real(double) :: zf(3,ns2)
    real(double) :: nzf(3,ns2,3,ns1,2)
    integer :: i,j,k,m

    real(double) :: theta,rot1(3,3),rot2(3,3),forceblk(3,3),mat1(3,3),mat2(3,3)
    integer :: ii,jj
    logical :: rotate_crystal

    write(*,'(A)') 'read: To open file ='//trim(abspartialfrcfile)
    open(unit=fu,file=trim(abspartialfrcfile),status='old',action='read')
    do i = 1, ns2
      read(fu,*) zf(1:3,i)
    enddo
    do i = 1, ns1
      do j = 1, 3
        do m = 1, 2
          do k = 1, ns2
            read(fu,*) nzf(1:3,k,j,i,m)
          enddo
        enddo
      enddo
    enddo
    close(fu)

    rotate_crystal = .false.

    if( .not. rotate_crystal) then

    else if(rotate_crystal) then

      theta = -1.0d0
      rot1(1:3,1:3) = reshape( (/ cos(theta), sin(theta),zero, -sin(theta),cos(theta),zero,zero,zero,one/), (/3,3/) )
      theta = -theta
      rot2(1:3,1:3) = reshape( (/ cos(theta), sin(theta),zero, -sin(theta),cos(theta),zero,zero,zero,one/), (/3,3/) )
      do m = 1, 2
        do i = 1, ns1
          do k = 1, ns2
            do ii = 1, 3
              do jj = 1 ,3
                forceblk(jj,ii) = nzf(jj,k,ii,i,m)
              enddo
            enddo

            call matmatn(3,rot1(1,1),forceblk(1,1),mat1(1,1))
            call matmatn(3,mat1(1,1),rot2(1,1),mat2(1,1))
            do ii = 1, 3
              do jj = 1, 3
                nzf(jj,k,ii,i,m) = mat2(jj,ii)
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
    write(*,'(A)') 'file='//trim(abspartialfrcfile)//' is now closed'
  end subroutine read_forces

  subroutine write_forces(ofu,absfinalforcesfile,ns1,ns2,zf,nzf)
    integer :: ofu
    character(len=*) :: absfinalforcesfile
    integer :: ns1,ns2
    real(double) :: zf(3,ns2)
    real(double) :: nzf(3,ns2,3,ns1,2)
    integer :: i,j,k,m

    write(*,'(A)') 'write: To open file ='//trim(absfinalforcesfile)
    open(unit=ofu,file=trim(absfinalforcesfile),status='replace')
    do i = 1, ns2
      write(ofu,*) zf(1:3,i)

    enddo
    do i = 1, ns1
      do j = 1, 3
        do m = 1, 2
          do k = 1, ns2

            write(ofu,*) nzf(1:3,k,j,i,m)

          enddo
        enddo
      enddo
    enddo
    close(ofu)
    write(*,'(A)') 'file='//trim(absfinalforcesfile)//' is now closed'
  end subroutine write_forces

  subroutine report_operation_type(dm_int,tr_int)
    integer :: dm_int,tr_int

    if(dm_int == 1 .and. tr_int == -1) then
      write(*,*) 'operation is 2 (2-fold rotation)'
    else if(dm_int == 1 .and. tr_int == 0) then
      write(*,*) 'operation is 3 (3-fold rotation)'
    else if(dm_int == 1 .and. tr_int == 1) then
      write(*,*) 'operation is 4 (4-fold rotation)'
    else if(dm_int == 1 .and. tr_int == 2) then
      write(*,*) 'operation is 6 (6-fold rotation)'
    else if(dm_int == 1 .and. tr_int == 3) then
      write(*,*) 'operation is 1 (identity)'

    else if(dm_int == -1 .and. tr_int == 1) then
      write(*,*) 'operation is -2 (mirror)'
    else if(dm_int == -1 .and. tr_int == 0) then
      write(*,*) 'operation is -3'
    else if(dm_int == -1 .and. tr_int == -1) then
      write(*,*) 'operation is -4'
    else if(dm_int == -1 .and. tr_int == -2) then
      write(*,*) 'operation is -6'
    else if(dm_int == -1 .and. tr_int == -3) then
      write(*,*) 'operation is -1 (inversion)'
    else
      write(*,*) 'determinant,trace= ',dm_int,tr_int
      write(*,'(A)') 'err: bad combination of trace and dmerminant.'
      stop 1
    endif
  end subroutine report_operation_type

  subroutine frac2nd(f,n,d)
    real(double) :: f,n,d,absf,intf,ffrac,p,intp,pfrac
    integer :: i
    integer :: nscan=1000
    real(double),parameter :: eps=1.0d-6
    logical :: found

    if(abs(f) < eps) then
      n = zero
      d = one
      return
    endif

    if(f > zero) then
      absf = f
    else
      absf = -f
    endif

    intf = nint(absf)
    ffrac = absf - intf

    if(abs(ffrac) < eps) then

      n = f
      d = one
      return
    endif

    found = .false.
    do i = 1, nscan
      p = i*ffrac

      intp = nint(p)
      pfrac = p - intp
      if(abs(pfrac) < eps) then
        n = intf*i + intp
        d = i
        if(f < zero) then
          n = -n
        endif

        if(abs(f - n/d) > eps) then
          write(*,*) 'nscan = ',nscan
          write(*,*) 'f,n,d,n/d,f-n/d=',f,n,d,n/d,f-n/d
          write(*,'(A)') 'err: wrong.'
          stop 1
        endif
        found = .true.
        return
      endif
    enddo
    if(.not. found) then
      write(*,*) 'we cannot find n and d for f=',f
      write(*,*) 'eps =',eps
      write(*,*) 'nscan = ',nscan
      write(*,'(A)') 'err: check the value of f or decrease eps or increase nscan?'
      stop 1
    endif
  end subroutine frac2nd

  subroutine LU(pivot,perm,n,A,L,U)
    logical :: pivot
    integer :: lprint
    integer :: n,i,j,k,tint,jmax,perm(n)
    real(double) :: absu,tmp,maxn,diff,m,A(n,n),L(n,n),U(n,n)
    real(double),allocatable :: T(:,:),PA(:,:)

    do i = 1, n
      perm(i) = i
    enddo

    do i = 1, n
      U(1:n,i) = A(1:n,i)
    enddo

    L(:,:) = zero

    do i = 1, n
      L(i,i) = one
      if(pivot) then

        maxn = 1.0d-100
        do j = i, n
          absu = abs(U(j,i))
          if(absu > maxn) then
            maxn = absu

            jmax = j
          endif
        enddo

        if(jmax /= i) then

          if(i < 1 .or. i > n) then
            write(*,*) 'i = ',i
            write(*,*) 'matrix may be singular?'
            write(*,'(A)') 'err: out of bound.'
            stop 1
          endif
          if(jmax < 1 .or. jmax > n) then
            write(*,*) 'jmax = ', jmax
            write(*,*) 'matrix may be singular?'
            write(*,'(A)') 'err: out of bound.'
            stop 1
          endif
          tint = Perm(i)
          Perm(i) = Perm(jmax)
          Perm(jmax) = tint

          do j = i, n
            tmp = U(i,j)
            U(i,j) = U(jmax,j)
            U(jmax,j) = tmp
          enddo

          do j = 1, i-1
            tmp = L(i,j)
            L(i,j) = L(jmax,j)
            L(jmax,j) = tmp
          enddo
        endif
      endif
      do j = i+1, n
        if(abs(U(i,i)) < 1.0d-10) then
          write(*,*) 'i,U(i,i) = ',i,U(i,i)
          write(*,'(A)') 'err: pivot value is zero.'
          stop 1
        endif
        m = U(j,i)/U(i,i)
        L(j,i) = m

        do k = i, n
          U(j,k) = U(j,k) - m*U(i,k)
        enddo
      enddo
    enddo

    lprint = 1
    if(lprint == 1) then
      write(*,*) 'Perm = '
      do i = 1, n
        write(*,'(I5)') Perm(i)
      enddo
      write(*,*) 'L = '
      do i = 1, n
        write(*,'(20F12.6)') L(i,1:n)
      enddo
      write(*,*) 'U = '
      do i = 1, n
        write(*,'(20F12.6)') U(i,1:n)
      enddo
      write(*,*)
    endif

    allocate(T(n,n))
    allocate(PA(n,n))

    call matmatn(n,L(1,1),U(1,1),T(1,1))

    do i = 1, n
      PA(i,1:n) = A(Perm(i),1:n)
    enddo

    diff = zero
    do i = 1, n
      do j = 1, n
        T(j,i) = T(j,i) - PA(j,i)
        diff = diff + abs(T(j,i))
      enddo
    enddo
    if(diff > 1d-10) then
      write(*,*) 'diff between PA and its LU  is ',diff
      write(*,'(A)') 'err: LU is not the same as PA.'
      stop 1
    else

    endif
    deallocate(T)
    deallocate(PA)
  end subroutine LU

  subroutine matrix_inverse(pivot,n,A,invA)
    logical :: pivot
    integer :: n
    real(double) :: A(n,n),invA(n,n)
    real(double),allocatable :: sol(:)
    real(double),allocatable :: L(:,:),U(:,:)
    real(double),allocatable :: newRHS(:),tempRHS(:)
    integer,allocatable :: Perm(:)
    integer :: i
    real(double) :: det

    allocate(L(n,n),U(n,n))
    allocate(Perm(n))
    allocate(newRHS(n))
    allocate(tempRHS(n))
    allocate(sol(n))

    call LU(pivot,Perm(1),n,A(1,1),L(1,1),U(1,1))

    det = one
    do i = 1, n
      det = det*U(i,i)
    enddo
    write(*,*) 'det(A) = ',det

    do i = 1, n
      tempRHS(1:n) = zero
      tempRHS(i) = one
      call solveX(pivot,n,L(1,1),U(1,1),Perm(1),tempRHS(1),newRHS(1),sol(1))
      invA(1:n,i) = sol(1:n)
    enddo

    write(*,*) 'inverse of A is '
    do i = 1, n
      write(*,'(100F15.5)') invA(i,1:n)
    enddo
    deallocate(L)
    deallocate(U)
    deallocate(Perm)
    deallocate(newRHS)
    deallocate(sol)
  end subroutine matrix_inverse

  subroutine solvelinear(pivot,n,A,RHS,sol)
    logical :: pivot
    integer :: n
    real(double) :: A(n,n)
    real(double) :: RHS(n)
    real(double) :: sol(n)
    real(double),allocatable :: L(:,:),U(:,:),invA(:,:)
    real(double),allocatable :: newRHS(:),tempRHS(:)
    integer,allocatable :: Perm(:)

    allocate(L(n,n),U(n,n))
    allocate(Perm(n))
    allocate(newRHS(n))
    allocate(tempRHS(n))
    allocate(invA(n,n))

    call LU(pivot,Perm(1),n,A(1,1),L(1,1),U(1,1))

    call solveX(pivot,n,L(1,1),U(1,1),Perm(1),RHS(1),newRHS(1),sol(1))

    if(pivot) then

    else

    endif

    deallocate(L)
    deallocate(U)
    deallocate(Perm)
    deallocate(newRHS)
    deallocate(invA)
  end subroutine solvelinear

  subroutine solveX(pivot,n,L,U,Perm,RHS,newRHS,sol)
    logical :: pivot
    integer :: n
    real(double) :: L(n,n),U(n,n)
    real(double),allocatable :: sol1(:),sol2(:)
    integer :: perm(n)
    real(double) :: RHS(n),newRHS(n),sol(n),s
    integer :: i,j
    integer,parameter :: no_print_opt=0
    integer,parameter :: print_opt=1
    integer :: print_int

    print_int=print_opt

    allocate(sol1(n),sol2(n))

    if(pivot) then

      do i = 1, n
        newRHS(i) = RHS(Perm(i))
      enddo
    else
      do i = 1, n
        newRHS(i) = RHS(i)
      enddo
    endif

    do i = 1, n
      s = zero
      do j = 1, i-1
        s = s + L(i,j)*sol1(j)
      enddo

      sol1(i) = newRHS(i) - s
    enddo

    if(print_int == print_opt) then
      write(*,'(A,10F10.5)') ' newRHS is ',newRHS(1:n)
      write(*,*) 'sol from L part is '
      write(*,'(20F20.6)') sol1(1:n)
      write(*,*)
    endif

    do i = n, 1, -1
      s = zero
      do j = n, i+1, -1
        s = s + U(i,j)*sol2(j)
      enddo
      sol2(i) = (sol1(i) - s)/U(i,i)
      sol(i) = sol2(i)
    enddo
    deallocate(sol1,sol2)
  end subroutine solveX

  subroutine coefxyz(coef,p1,p2,BTB)
    real(double) :: coef(4),p1(3),p2(3),BTB(3,3)
    real(double) :: p,q,r, t,u,v, C11,C12,C13,C21,C22,C23,C31,C32,C33

    p = p1(1); q = p1(2); r = p1(3)
    t = p2(1); u = p2(2); v = p2(3)

    C11 = BTB(1,1); C12 = BTB(1,2); C13 = BTB(1,3)
    C21 = BTB(2,1); C22 = BTB(2,2); C23 = BTB(2,3)
    C31 = BTB(3,1); C32 = BTB(3,2); C33 = BTB(3,3)

    coef(1) = (t-p)*2*C11 + (u-q)*C12 + (v-r)*C13 + (u-q)*C21 + (v-r)*C31
    coef(2) = (t-p)*C12 + (t-p)*C21 + (u-q)*2*C22 + (v-r)*C23 + (v-r)*C32
    coef(3) = (t-p)*C13 + (u-q)*C23 + (t-p)*C31 + (u-q)*C32 + (v-r)*2*C33
    coef(4) = (t-p)*(-p-t)*C11 + (p*q - t*u)*C12 + (p*r - t*v)*C13 + (q*p - u*t)*C21 + (u-q)*(-q-u)*C22 + (q*r - u*v)*C23 + (r*p - v*t)*C31 + (r*q - v*u)*C32 + (v-r)*(-r-v)*C33

  end subroutine coefxyz

  subroutine kptxyz(BTB,sh1,sh2,th1,th2,uh1,uh2,sol)
    real(double) :: BTB(3,3)
    real(double) :: sh1(3),sh2(3)
    real(double) :: th1(3),th2(3)
    real(double) :: uh1(3),uh2(3)
    real(double) :: mat(3,3),b(3),sol(3),coef(4)
    integer :: iprint

    call coefxyz(coef(1),sh1(1),sh2(1),BTB(1,1))
    mat(1,1:3) = coef(1:3)
    b(1) = -coef(4)

    call coefxyz(coef(1),th1(1),th2(1),BTB(1,1))
    mat(2,1:3) = coef(1:3)
    b(2) = -coef(4)

    call coefxyz(coef(1),uh1(1),uh2(1),BTB(1,1))
    mat(3,1:3) = coef(1:3)
    b(3) = -coef(4)

    iprint = 1
    if(iprint == 1) then
      write(*,*) 'mat = '
      write(*,'(3F20.10)') mat(1,1:3)
      write(*,'(3F20.10)') mat(2,1:3)
      write(*,'(3F20.10)') mat(3,1:3)
      write(*,*) 'b = '
      write(*,'(3F20.10)') b(1:3)
    endif

    call solvelinear(.true.,3,mat(1,1),b(1),sol(1))

    if(iprint == 1) then
      write(*,*) 'sol is '
      write(*,'(20F12.6)') sol(1:3)
    endif
  end subroutine kptxyz

  subroutine form_groups(nk,rsq,ig,ng,degentol)
    integer :: nk
    real(double) :: rsq(nk)
    integer :: ig(nk+1),ng,i,subdim
    real(double) :: degentol
    real(double) :: rg
    integer :: tsum

    rg = rsq(1)
    ng = 1
    ig(1) = 1
    do i = 2, nk
      if(abs(rsq(i)-rg) > degentol) then
        rg = rsq(i)
        ng = ng + 1
        ig(ng) = i
      endif
    enddo
    ig(ng+1) = nk+1

    tsum = 0

    do i = 1, ng
      subdim = ig(i+1)-ig(i)
      tsum = tsum + subdim

    enddo
    if(tsum /= nk) then
      write(*,*) 'tsum,nk=',tsum,nk
    endif
  end subroutine form_groups

  subroutine pt_shiftedfreq(aev,pev,zfreq,pfreq,nstages,minstage,degentol,nb,D1,D2,deltaD,tmpv1,tmpv2,tmpv3,T1,T2,w1,w2,indarr,CWORK,LCWORK,RWORK,LRWORK,Di1,Di2)

    integer :: nb,nstages,stage,LCWORK,LRWORK,minstage
    integer :: ng,m,k,dimen,mm,info,indarr(nb)
    integer :: dp,ind,ind2,kk,indw2,recordkk,tind

    integer,allocatable :: localmap1(:),localmap2(:),ig(:),gmap(:,:)
    integer,allocatable :: ipiv(:),store(:)

    integer,parameter :: ndata = 3
    integer,parameter :: ndeg=1

    complex(double) :: CWORK(LCWORK),D1(nb,nb),D2(nb,nb),deltaD(nb,nb)
    complex(double) :: tmpv1(nb,nb),tmpv2(nb,nb)
    complex(double) :: tmpv3(nb,nb)
    complex(double) :: T1(nb,nb),T2(nb,nb)
    complex(double) :: Di1(nb,nb),Di2(nb,nb)
    complex(double) :: Calpha,Cbeta

    real(double) :: pev(nb,0:nstages),aev(nb,0:nstages)
    real(double) :: offset,newv,mindis,dis,newvp,freqsq
    real(double) :: RWORK(LRWORK),zfreq(nb),pfreq(nb),w1(nb),w2(nb),freq,degentol

    real(double),allocatable :: tmpnu(:),deltawsqr(:),predictE(:),lambda(:),xval(:),yval(:)
    real(double),allocatable :: Mv(:,:),MT(:,:),MTM(:,:),MTy(:),sol(:,:),s(:)
    real(double),allocatable :: freqsqarr(:,:),sortarr(:),sev(:,:)
    integer :: debug

    allocate(xval(ndata))
    allocate(yval(ndata))
    allocate(MV(ndata,0:ndeg))
    allocate(MT(0:ndeg,ndata))
    allocate(MTM(0:ndeg,0:ndeg))
    allocate(MTy(0:ndeg))
    allocate(sol(0:ndeg,1))
    allocate(s(0:ndeg))
    allocate(ipiv(0:ndeg))

    allocate(store(nb))

    allocate(lambda(0:nstages))
    allocate(freqsqarr(nb,0:nstages))

    allocate(gmap(nb,0:nstages))
    allocate(ig(nb+1))
    allocate(deltawsqr(nb))
    allocate(predictE(nb))
    allocate(localmap1(nb))
    allocate(localmap2(nb))

    allocate(sortarr(nb))
    allocate(sev(nb,0:nstages))
    allocate(tmpnu(nb))

    if(nstages < 1) then
      write(*,*) 'nstages =',nstages
      stop 1
    endif

    do stage = 0, nstages
      lambda(stage) = stage*1.0D0/(nstages*1.0D0)
    enddo

    tmpv1(1:nb,1:nb) = D1(1:nb,1:nb)
    call zheev('V','U',nb,tmpv1(1,1),nb,w1(1),CWORK(1),LCWORK,RWORK(1),info)
    if(info /= 0) then
      write(*,*) 'info = ',info
      write(*,'(A)') 'err: zheev error.'
      stop 1
    endif
    freqsqarr(1:nb,0) = w1(1:nb)

    do k = 1, nb
      gmap(k,0) = k
    enddo

    aev(1:nb,0) = w1(1:nb)
    pev(1:nb,0) = w1(1:nb)

    do k = 1, nb
      freqsq = w1(k)
      freq = sqrt(abs(freqsq))
      if(freqsq < zero) then
        tmpnu(k) = -freq*Rydoverh2invcm
      else
        tmpnu(k) = freq*Rydoverh2invcm
      endif
    enddo

    debug = 0
    if(debug == 1) then
      write(*,'(A)') 'eigenvalues (in 1/cm) of zeroth structure, w1, are:'
      call print_arr_in_n_columns(tmpnu(1),nb,6)
    endif

    w2(1:nb) = w1(1:nb)
    tmpv2(1:nb,1:nb) = tmpv1(1:nb,1:nb)
    Di2(1:nb,1:nb) = D1(1:nb,1:nb)

    do k = 1, nb
      do m = 1, nb
        deltaD(m,k) = (D2(m,k) - D1(m,k))/(nstages*one)
      enddo
    enddo

    do stage = 1, nstages

      w1(1:nb) = w2(1:nb)
      tmpv1(1:nb,1:nb) = tmpv2(1:nb,1:nb)
      Di1(1:nb,1:nb) = Di2(1:nb,1:nb)

      Di2(1:nb,1:nb) = D1(1:nb,1:nb) + dcmplx(lambda(stage),0)*(D2(1:nb,1:nb) - D1(1:nb,1:nb))
      tmpv2(1:nb,1:nb) = Di2(1:nb,1:nb)

      call zheev('V','U',nb,tmpv2(1,1),nb,w2(1),CWORK(1),LCWORK,RWORK(1),info)
      if(info /= 0) then
        write(*,*) 'info = ',info
        write(*,'(A)') 'err: zheev error.'
        stop 1
      endif
      freqsqarr(1:nb,stage) = w2(1:nb)

      predictE(1:nb) = zero

      if(stage < minstage) then

        call form_groups(nb,w1(1),ig(1),ng,degentol)

        do k = 1, ng
          dimen = ig(k+1) - ig(k)

          do m = 1, dimen
            do mm = 1, nb
              tmpv3(mm,m) = tmpv1(mm,ig(k)+(m-1))
            enddo
          enddo

          T1(1:nb,1:dimen) = CZero
          Calpha = COne
          Cbeta = CZero
          call zgemm('N','N',nb,dimen,nb,Calpha,deltaD(1,1),nb,tmpv3(1,1),nb,Cbeta,T1(1,1),nb)
          T2(1:dimen,1:dimen) = CZero
          Calpha = COne
          Cbeta = CZero
          call zgemm('C','N',dimen,dimen,nb,Calpha,tmpv3(1,1),nb,T1(1,1),nb,Cbeta,T2(1,1),nb)

          call zheev('N','U',dimen,T2(1,1),nb,deltawsqr(1),CWORK(1),LCWORK,RWORK(1),info)
          if(info /= 0) then
            write(*,*) 'info = ',info
            write(*,'(A)') 'err: zheev error.'
            stop 1
          endif

          do m = 1, dimen
            predictE(ig(k)+(m-1)) = w1(ig(k)+(m-1)) + deltawsqr(m)
          enddo
        enddo

        do k = 1, nb
          indarr(k) = k
        enddo
        call dbl_sort(nb,predictE(1),indarr(1),1)

        call inversemapping(nb,indarr(1),localmap1(1))

      else

        do k = 1, nb
          store(k) = k
        enddo

        do k = 1, nb

          do dp = 1, ndata
            xval(dp) = lambda(stage-ndata+dp-1)
            yval(dp) = aev(k,stage-ndata+dp-1)
          enddo
          offset = lambda(stage)
          call least_sq_fit(ndeg,ndata,offset,xval(1),yval(1),s(0),Mv(1,0),MT(0,1),MTM(0,0),MTy(0),sol(0,1),ipiv(0))

          call fx_fpx_least_sq_fit(lambda(stage),s(0),ndeg,offset,newv,newvp)

          mindis = 1.0d100
          recordkk = -1
          do kk = nb, k, -1
            indw2 = store(kk)
            dis = abs(w2(indw2)-newv)
            if(dis < mindis) then

              mindis = dis
              recordkk = kk
            endif
          enddo

          if(recordkk == -1) then
            write(*,*) 'Bad recordkk= ',recordkk
            stop 1
          endif
          localmap2(k) = store(recordkk)

          tind = store(k)
          store(k) = store(recordkk)
          store(recordkk) = tind

        enddo

      endif

      do k = 1, nb
        ind = gmap(k,stage-1)
        if(stage < minstage) then
          ind2 = localmap1(ind)
          gmap(k,stage) = ind2
        elseif(stage >= minstage) then
          gmap(k,stage) = localmap2(k)
        endif
      enddo

      do k = 1, nb
        pev(k,stage) = predictE(gmap(k,stage))
        aev(k,stage) = w2(gmap(k,stage))
      enddo
    enddo

    do stage = nstages-ndata, 0, -1

      do k = 1, nb
        store(k) = k
      enddo

      do k = 1, nb
        do dp = 1, ndata
          xval(dp) = lambda(stage+dp)
          yval(dp) = aev(k,stage+dp)
        enddo
        offset = lambda(stage)
        call least_sq_fit(ndeg,ndata,offset,xval(1),yval(1),s(0),Mv(1,0),MT(0,1),MTM(0,0),MTy(0),sol(0,1),ipiv(0))
        call fx_fpx_least_sq_fit(lambda(stage),s(0),ndeg,offset,newv,newvp)
        mindis = 1.0d100

        do kk = nb, k, -1
          indw2 = store(kk)
          dis = abs(freqsqarr(indw2,stage)- newv)
          if(dis < mindis) then
            mindis = dis
            recordkk = kk
          endif
        enddo

        aev(k,stage) = freqsqarr(  store(recordkk), stage )

        tind = store(k)
        store(k) = store(recordkk)
        store(recordkk) = tind

      enddo
    enddo

    do k = 1, nb
      indarr(k) = k
      sortarr(k) = aev(k,0)
    enddo
    call dbl_sort(nb,sortarr(1),indarr(1),1)
    do k = 1, nb
      do stage = 0, nstages
        sev(k,stage) = aev(indarr(k),stage)
      enddo
    enddo

    do k = 1, nb
      do stage = 0, nstages
        aev(k,stage) = sev(k,stage)
      enddo
    enddo

    do k = 1, nb
      freqsq = aev(k,0)
      freq = sqrt(abs(freqsq))
      if(freqsq < zero) then
        zfreq(k) = -freq*Rydoverh2invcm
      else
        zfreq(k) = freq*Rydoverh2invcm
      endif
    enddo

    do k = 1, nb
      freqsq = aev(k,nstages)
      freq = sqrt(abs(freqsq))
      if(freqsq < zero) then
        pfreq(k) = -freq*Rydoverh2invcm
      else
        pfreq(k) = freq*Rydoverh2invcm
      endif
    enddo

    deallocate(xval)
    deallocate(yval)
    deallocate(Mv)
    deallocate(MT)
    deallocate(MTM)
    deallocate(MTy)
    deallocate(s)
    deallocate(sol)
    deallocate(ipiv)

    deallocate(store)

    deallocate(lambda)
    deallocate(freqsqarr)

    deallocate(gmap)
    deallocate(ig)
    deallocate(deltawsqr)
    deallocate(predictE)
    deallocate(localmap1)
    deallocate(localmap2)

    deallocate(sortarr)
    deallocate(sev)
    deallocate(tmpnu)

  end subroutine pt_shiftedfreq

  subroutine find_fit_svd(m,n,a,w,v,y,c,f)
    integer :: m,n
    real(double) :: a(m,n),w(n),v(n,n),y(m),c(n),f(n)
    integer :: i,j

    do i = 1, n
      c(i) = 0
      do j = 1, m
        c(i) = c(i) + a(j,i)*y(j)
      enddo
      write(*,*) 'i,c(i)=', i, c(i)
    enddo

    write(*,*) 'now z:'
    do i = 1, n
      c(i) = c(i)/w(i)
      write(*,*) 'i,z(i)=', i, c(i)
    enddo

    write(*,*)

    do i = 1, n
      f(i) = 0
      do j = 1, n
        f(i) = f(i) + v(i,j)*c(j)
      enddo
      write(*,*) 'i,f(i)=', i, f(i)
    enddo
  end subroutine find_fit_svd

  subroutine bco_para_from_dht_strained_trigonal(d,h,t)
    real(double) :: g,d,h,t,a,b,c,cosg,cosg2,ba,ca,bc,bs

    a = 3d0*h
    b = sqrt(4*h*h+t*t)
    c = 2d0*d
    cosg = -2d0*h/b
    ba = b/a
    ca = c/a
    g = atan2(t,-2d0*h)
    cosg2 = cos(g)
    if(abs(cosg2-cosg) > 1d-10) then
      write(*,'(A)') 'err: cos diff.'
      stop 1
    endif
    bc = b*cosg
    bs = b*sin(g)

    write(*,*) 'd = ',d
    write(*,*) 'h = ',h
    write(*,*) 't = ',t

    write(*,*)
    write(*,*) 'a = ',a
    write(*,*) 'For potential QE use: a*Ang2Bohr = ',a*Ang2Bohr

    write(*,*) 'b/a = ',ba
    write(*,*) 'c/a = ',ca
    write(*,*) 'cos(g) = ',cosg
    write(*,*)

    write(*,*) 'A1: ibrav0 (row view):'
    write(*,'(3F18.12)') d, h/2,t
    write(*,'(3F18.12)') -d,h/2,t

    write(*,'(3F18.12)') zero,-h,t

    write(*,*) 'A4: Before centering: QE orthorhombic form (row view):'
    write(*,'(3F18.12)') 3*h,zero,zero
    write(*,'(3F18.12)') -2*h, t, zero
    write(*,'(3F18.12)') zero, zero, 2*d

    write(*,*) 'A5: base-center monoclinic form: (row view)'
    write(*,'(3F18.12)') 3*h/2,zero,-d
    write(*,'(3F18.12)') -2*h, t, zero
    write(*,'(3F18.12)') 3*h/2, zero, d
  end subroutine bco_para_from_dht_strained_trigonal

  function gcd_i(a,b)
    integer :: gcd_i,a,b
    integer :: tempn,n1,n2,m,r
    logical :: done
    real(double) :: mreal

    n1 = abs(a)
    n2 = abs(b)

    if(n2 < n1) then
      tempn = n1
      n1 = n2
      n2 = tempn
    endif

    if(n1 == 0) then
      write(*,'(A)') 'err: n1 = 0.'
      stop 1
    endif

    done = .false.
    do
      if(done) exit

      r = modulo(n2,n1)

      mreal = (n2*one - r*one)/(n1*one)
      m = nint(mreal)
      if(r == 0) then
        done = .true.
        gcd_i = n1
      else

        n2 = n1
        n1 = r
      endif
    enddo
  end function gcd_i

  function lcm_i(a,b)
    integer :: lcm_i,a,b
    integer :: n1,n2,gcd
    real(double) :: reallcm

    n1 = abs(a)
    n2 = abs(b)

    gcd = gcd_i(n1,n2)
    reallcm = n1*one*n2/(gcd*one)
    lcm_i = nint(reallcm)
  end function lcm_i

  function lcm_arr(a,n)
    integer :: lcm_arr,n,a(n)
    integer :: i,lcmnow,newlcm
    if(n < 2) then
      write(*,'(A)') 'err: n is too small.'
      stop 1
    endif
    lcmnow = lcm_i(a(1),a(2))
    do i = 3, n
      newlcm = lcm_i(lcmnow,a(i))
      lcmnow = newlcm
    enddo
    lcm_arr = lcmnow
  end function lcm_arr

  subroutine print_arr_in_n_columns(r,nb,n)
    integer :: nb,n
    real(double) :: r(nb)
    integer :: k,rem

    if(n < 1) then
      write(*,*) 'In print_arr_in_n_columns: n=',n
      write(*,*) 'n must be greater than 1'
      stop 1
    endif
    do k = 1, nb/n
      write(*,'(50E18.8)') r( (k-1)*n+1: n*k)
    enddo

    rem = modulo(nb,n)
    if(rem > 0) then
      write(*,'(50E18.8)') r( (nb/n)*n+1: (nb/n)*n+rem)
    endif
  end subroutine print_arr_in_n_columns

  subroutine reformat_file(filename)
    character(len=*) :: filename
    character(len=DSL) :: tmpstr
    integer :: ln,i
    integer :: ou
    integer :: fu

    fu = newunit()
    ln = num_of_lines(fu,filename)
    write(*,*) 'ln = ',ln

    fu = newunit()
    open(unit=fu,file=trim(filename),status='old',action='read')

    ou = newunit()
    open(unit=ou,file='out',status='replace')
    do i = 1, ln
      read(fu,DSF) tmpstr
      write(ou,'(A)') trim(tmpstr)
    enddo
    close(fu)
  end subroutine reformat_file

end module commod
