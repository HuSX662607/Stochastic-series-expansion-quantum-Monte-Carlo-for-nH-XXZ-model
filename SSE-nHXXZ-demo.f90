! Demo code of SSE-QMC for nH-XXZ model                                            !
! This is the source code for the article 'Nontrivial worldline winding in         !
! non-Hermitian quantum systems'(arXiv:2307.01260). The code is partially based on 
! the code written by Sandvik (https://physics.bu.edu/~sandvik/programs/ssebasic/ssebasic.f90)
! and the realization of some update methods are inspired by GitHub open-source    !
! codes https://github.com/joaogci/SSE/tree/master.
!                                                                                  !
! Labeling of vertices:                                                            !
!                                                                                  !
!  ↑     ↓       ↓     ↑       ↑     ↓       ↓     ↑       ↑     ↑       ↓     ↓   !
! |=======|     |=======|     |||||||||     |||||||||     |=======|     |=======|  !
!  ↑     ↓       ↓     ↑       ↓     ↑       ↑     ↓       ↑     ↑       ↓     ↓   !
!                                                                                  !
!    a=1           a=2           a=3           a=4           a=5           a=6     !


!--------------------!
 module configuration
!--------------------!
 save

 integer :: nn       ! number of sites;
 integer :: nb       ! number of bonds; nb=nn for PBC
 integer :: if_PBC   ! 1=PBC,0=OBC

 integer,parameter :: msteps=20000
 integer,parameter :: nbins=3
 integer,parameter :: isteps=40000

 integer :: nh       ! number of nontrivial operators
 integer :: mm       ! cutoff
 integer :: Nloop    ! number of loops
 integer :: W_Loop=0   ! worldline winding number
 integer :: total_loopsize=0
 integer :: N_WSector(0:nbins-1,-40:40)
 integer :: N_spinsector(0:nbins-1,-60:60)

 integer, allocatable :: bsites(:,:)   ! Array representing lattice geometry, i.e. how the sites are connected.
                                       ! index1: 0 or 1; index2: bond label b
                                       ! output: site label i

 integer, allocatable :: frstspinop(:) ! first operation on each site in linked vertex list
                                       ! index: site label i
                                       ! output: first nonnegative leg

 integer, allocatable :: lastspinop(:) ! last operation on each site in linked vertex list
                                       ! index: site label i
                                       ! output: last nonnegative leg

 integer, allocatable :: vertexlist(:) ! Array recording how the vertices are connected in a configuration
                                       ! index: leg label
                                       ! output: leg that connects the index leg

 integer, allocatable :: spin(:)       ! Array of spin configuration
                                       ! index: site label
                                       ! output: spin state (rep. by -1/1)

 integer, allocatable :: opstring(:)   ! Array of operator string
                                       ! index: time slice label p
                                       ! output: operator label 6*(b+1)+a-1

 real(8),dimension(0:5,0:3,0:3) :: VertexInfo ! Array recording information of vertices
                                              ! index1: vertex label a-1, a from 1 to 6; index2: entrance leg; index3: exit leg
                                              ! output: new vertex label a_new after flipping the spins of in&out legs

 real(8),dimension(0:5) :: MatrixElement   ! Matrix elements of vertices
                                           ! index: vertex label a-1
                                           ! output: matrix elements of the corresponding vertices

 real(8) :: Jz
 real(8) :: dJ
 real(8) :: delta
 real(8) :: beta
 real(8) :: mu
 real(8) :: C=0.5


 end module configuration
!------------------------!


!==============!
 program SSERun
!==============!
 use configuration;implicit none
 include 'mpif.h'
 integer :: i,j,m,k,wp,ierr,my_rank,p
 integer :: npc,ntotal

 real(8) :: mu_out,delta_out,E
 real(8) :: S,EL(0:nbins-1),EAndVar(2)
 character*8 :: str

write(str,'(I4)') my_rank

open(unit=10,file='Results_PBC_varJz'//trim(adjustl(str))//'.txt',status='replace')


mu_out=0
delta_out=0.2

call main(mu_out,delta_out,EL)

write(10,*) 'Jz=',Jz,'  dJ=',dJ,'  delta=',delta,'  -E/N=',EAndVar(1),'  variation:',EAndVar(2)
write(10,*) 'wp=',wp
write(10,*) N_WSector


close(10)
 end program SSERun
!=====================!


!----------------!
subroutine main(mu_out,delta_out,EAndVar)
!----------------!
use configuration;implicit none
integer :: i,j,k,wp,wpL(1)
real(8) :: mu_out,delta_out
real(8) :: E,Eavg,EAndVar(2),EList(nbins)

N_WSector(:,:)=0
N_spinsector(:,:)=0
nn=64
Jz=0
delta=delta_out
dJ=0.3
mu=mu_out
beta=100
! delta=0.5
if_PBC=1
 C=0.75
MatrixElement=(/ C+Jz*0.25d0, 0.5d0*(1-delta), C+Jz*0.25d0, 0.5d0*(1+delta), C-Jz*0.25d0+mu*0.5, C-Jz*0.25d0-mu*0.5/)
 Nloop=nn/4
Eavg=0d0

call FillVertex()
call makelattice1D()
call initconfig()


do i=1,msteps
    call diagonalupdate()
    call LinkedVertex()
    call LoopUpdate()
    call adjustcutoff(i)
    ! write(*,*) i,'mm=',mm
enddo
do k=0,nbins-1
    E=0
    do j=1,isteps
        call diagonalupdate()
        call LinkedVertex()
        if (if_PBC==1) then
            call LoopUpdate()
        elseif (if_PBC==0) then
            call LoopUpdate()
        endif
        call MeasureWinding()
        ! write(*,*) j,'w=',W_Loop
        E=E+dble(nh)/dble(beta)
        N_WSector(k,W_Loop)=N_WSector(k,W_Loop)+1
        if (abs(sum(spin))<=60) then
            N_spinsector(k,sum(spin)) = N_spinsector(k,sum(spin))+1
        endif
    enddo
    E=E/dble(nn*j)-dble(C*nb)/dble(nn)
    EList(k+1)=E
    Eavg=Eavg+E
enddo

Eavg=Eavg/3.0d0

EAndVar(1)=Eavg
EAndVar(2)=sqrt(((EList(1)-Eavg)**2+(EList(2)-Eavg)**2+(EList(3)-Eavg)**2)/3.0d0)

wpL=maxloc(N_WSector(nbins-1,:))
wp=wpL(1)-41

write(*,*) 'Jz=',Jz,'  dJ=',dJ,'  delta=',delta,'μ=',mu,'  -E/N=',EAndVar(1),'  variation:',EAndVar(2)
write(*,*) 'wp=',wp
write(*,*) 'w:'
write(*,*) N_WSector(0,:)
write(*,*) N_WSector(1,:)
write(*,*) N_WSector(2,:)
write(*,*) 'spin distribution:'
write(*,*) N_spinsector(0,:)
write(*,*) N_spinsector(1,:)
write(*,*) N_spinsector(2,:)

call deallocateall()

return

end subroutine main
!------------------!

!----------------------!
subroutine FillVertex()                      
!----------------------!
! Filling the information of vertices into the array Vertexinfo !
! vertex leg label:                                             !
!                       2     3                                 !
!                      |=======|                                !
!                       0     1                                 !
use configuration;implicit none
integer,dimension(0:3) :: vertexspin,vertexspin_aft
integer :: a,xin,xout

do a=1,6
    a=a-1
    if ((a==0) .or. (a==3) .or. (a==4)) then
    ! Check the current spin states of the vertex a
        vertexspin(0) = 1
        if (a==4) then
            vertexspin(1) = 1
        else
            vertexspin(1) = -1
        endif
    else
        vertexspin(0) = -1
        if (a==5) then
            vertexspin(1) = -1
        else
            vertexspin(1) = 1
        endif
    endif
    if ((a==0) .or. (a==2) .or. (a==4) .or. (a==5)) then
        vertexspin(2) = vertexspin(0)
        vertexspin(3) = vertexspin(1)
    else
        vertexspin(2) = -vertexspin(0)
        vertexspin(3) = -vertexspin(1)
    endif
    vertexspin_aft(:)=vertexspin(:)
    do xin=0,3
        do xout=0,3
        ! Recording the resulting vertex index after flipping the spins at xin (entrance leg label) and xout (exit leg label)
            vertexspin_aft(xin) = -vertexspin_aft(xin)
            vertexspin_aft(xout) = -vertexspin_aft(xout)
            if (vertexspin_aft(0)*vertexspin_aft(1)>0) then
                if (vertexspin_aft(0) == vertexspin_aft(2)) then
                    VertexInfo(a,xin,xout) = int((11-vertexspin_aft(0))/2)
                else
                    VertexInfo(a,xin,xout) = -1
                endif
            elseif (vertexspin_aft(0) == vertexspin_aft(2)) then
                VertexInfo(a,xin,xout) = int(2-vertexspin_aft(0))
            else
                VertexInfo(a,xin,xout) = int(3+vertexspin_aft(0))
            endif
            vertexspin_aft(:)=vertexspin(:)
        enddo
    enddo
enddo
end subroutine FillVertex
!----------------------!


!----------------------!
subroutine initconfig()                  
!----------------------!
! Subroutine generating initial configuration !
use configuration;implicit none
integer :: i,j,i_A,i_A_large,i_A_small,i_A_d,ip,x
allocate(spin(0:nn-1))

mm=10

 allocate(opstring(0:mm-1))
 allocate(frstspinop(0:nn-1))
 allocate(lastspinop(0:nn-1))
 allocate(vertexlist(0:4*mm-1))
 opstring(:)=0
 vertexlist(:)=0
 nh=0

 spin(:)=-1
 do i=0,nn/2
    x=int(ran()*nn)
    do while (spin(x)==1)
        x=int(ran()*nn)
    enddo
    spin(x)=1
 enddo

end subroutine initconfig
!------------------------!


!--------------------------!
 subroutine makelattice1D()           
!--------------------------! 
! Subroutine generating the lattice geometry !
 use configuration; implicit none

 integer :: s,x1,x2,y1,y2

 nb=nn-1+if_PBC

 allocate(bsites(0:1,0:nb-1))

 
 do x1=0,nn-2
    bsites(0,x1)=x1
    bsites(1,x1)=x1+1       
 enddo

 if (if_PBC==1) then
    bsites(0,nn-1)=nn-1
    bsites(1,nn-1)=0
 endif

 end subroutine makelattice1D
!----------------------------!


!-------------------------!
subroutine diagonalupdate()
!-------------------------!
! Subroutine of diagonal updates !

use configuration;implicit none
integer,dimension(0:nn-1) :: spin0_1
integer :: p,i_op,op,b,a,i,ip,j_A
real(8) :: rd,Pa,Pd,aprob,dprob

spin0_1(:)=spin(:)
do p=0,mm-1
    op = opstring(p)
    if (op==0) then
    ! No operator at time slice p: insert a diagonal one
        b = int(ran()*nb)
        rd = ran()
        Pa = rd*(mm-nh)
        if (Pa < nb*beta*(C-spin(bsites(0,b))*spin(bsites(1,b))*Jz*0.25d0+mu*0.5*0.5*(spin(bsites(0,b))+spin(bsites(1,b))))) then
            if (spin(int(bsites(0,b)))==1) then
                if (spin(int(bsites(1,b)))==-1) then
                ! spins antiparallel ↑↓
                    a = 1
                elseif (spin(int(bsites(1,b)))==1) then
                ! spins parallel ↑↑
                    a = 5                        
                endif
            elseif (spin(int(bsites(0,b)))==-1) then
                if (spin(int(bsites(1,b)))==-1) then
                ! spins parallel ↓↓
                    a = 6
                elseif (spin(int(bsites(1,b)))==1) then
                ! spins antiparallel ↓↑
                    a = 3
                endif
            endif
            opstring(p)=6*(b+1) + a - 1
            nh = nh+1
            if ((spin(int(bsites(0,b)))*spin(int(bsites(1,b)))>0) .and. ((a==1) .or. (a==3))) then
            ! check if the configuration is valid for the model
                write(*,*) 'invalid operator'
            elseif ((spin(int(bsites(0,b)))*spin(int(bsites(1,b)))<0) .and. ((a==5) .or. (a==6))) then
                write(*,*) 'invalid operator'
            endif
        endif

    elseif (mod(op,2)==0 .or. (int(mod(op,6)+1)==6)) then
    ! diagonal operator at slice p: deleting the operator
        b = int(op/6 -1)
        a = int(mod(op,6)+1)
        rd = ran()
        Pd = rd/dble(mm-nh+1)
        aprob = nb*beta*(C-spin(bsites(0,b))*spin(bsites(1,b))*Jz*0.25d0+mu*0.5*0.5*(spin(bsites(0,b))+spin(bsites(1,b))))
        dprob = 1/aprob
        if (Pd<dprob) then
            opstring(p) = 0
            nh = nh-1
        endif
        if ((spin(int(bsites(0,b)))*spin(int(bsites(1,b)))>0) .and. ((a==1) .or. (a==3))) then
            write(*,*) 'invalid operator'
        elseif ((spin(int(bsites(0,b)))*spin(int(bsites(1,b)))<0) .and. ((a==5) .or. (a==6))) then
            write(*,*) 'invalid operator'
        endif

    else
    ! off-diagonal operator at slice p: propagating the spin state
        b = int(op/6 -1)
        spin(int(bsites(0,b))) = -spin(int(bsites(0,b)))
        spin(int(bsites(1,b))) = -spin(int(bsites(1,b)))
        if ((spin(int(bsites(0,b)))*spin(int(bsites(1,b)))>0)) then
            write(*,*) 'invalid operator'
        endif
    endif
enddo
do i=0,nn-1
    if (spin0_1(i) /= spin(i)) then
    ! check if the configuration is valid for the model
        write(*,*) 'time period is broken'
        exit
    endif
enddo
end subroutine diagonalupdate
!----------------------------!


!------------------------!
subroutine LinkedVertex()
!------------------------!
! Generating the linked vertex representation. Leg v=4p+l where l is the !
! leg label for a single vertex (0~3). If the operator at p is identity, !
! v=-1.                                                                  !
use configuration;implicit none
integer :: v0,v1,v2,op,b,s1,s2,ip,j_A

frstspinop(:) = -1
lastspinop(:) = -1

do v0=0,4*mm-1,4
    op=opstring(v0/4)
    if (op/=0) then
        b = op/6 - 1
        s1 = bsites(0,b)
        s2 = bsites(1,b)
        v1 = lastspinop(s1)
        v2 = lastspinop(s2)
        if (v1 /= -1) then
            vertexlist(v1) = v0
            vertexlist(v0) = v1
        else
            frstspinop(s1) = v0
        endif
        if (v2 /= -1) then
            vertexlist(v2) = v0+1
            vertexlist(v0+1) = v2
        else
            frstspinop(s2) = v0+1
        endif
        lastspinop(s1) = v0+2
        lastspinop(s2) = v0+3
    else
        vertexlist(v0:v0+3) = -1
    endif
enddo 
do s1=0,nn-1
    v1 = frstspinop(s1)
    if (v1 /= -1) then
        v2 = lastspinop(s1)
        vertexlist(v1)=v2
        vertexlist(v2)=v1
    endif
enddo

end subroutine LinkedVertex
!--------------------------!


!---------------------!
subroutine LoopUpdate()
!---------------------!
use configuration;implicit none
integer :: a,b,innerv,v0,v1,v2,k,p,ip,i,j,asubs,j_A,int_time1,int_time2,loop_size
integer :: IfMinus,IfOff,stringcopy(0:mm-1)
real(8) :: r,W(0:3)


k=0
do
    if (nh==0) then
        exit
    endif
    v0 = int(rand()*(4*mm)) 
    ! Randomly pick a initial leg for directed loop update

    stringcopy(:)=opstring(:)
    call system_clock(int_time1)
    loop_size=0
    do
        if (opstring(v0/4)==0) then
        ! The leg can't belong to identity operator
            v0 = int(rand()*(4*mm))
        else
            exit
        endif
    enddo
    if (vertexlist(v0)<0) then
        cycle
    endif
    v1=v0

    do
        innerv=int(mod(v1,4))
        ! the leg label (0~3) for the v1, as the entrance leg of the vertex

        p=v1/4
        b=opstring(p)/6-1
        a=int(mod(opstring(p),6))+1
        r = rand()
        if (r==0) then
            r=rand()
        endif
        W(:)=0
        do j=0,3
            W(j) = W(j)+IfMinus(a,innerv,j)*(MatrixElement(int(VertexInfo(a-1,innerv,j)-1)))
            W(j) = W(j)+IfMinus(a,innerv,j)*(- IfOff(a,innerv,j)*0.5d0*dJ * dble((-1)**bsites(0,b)))
        enddo
        if ((r>0) .and. (r<=W(0)/sum(W))) then
        ! the exit leg is chosen as 0 by Monte-Carlo
            if (W(0)==0) then
                write(*,*) 'wrong element'
            endif
            v2 = 4*p + 0
            asubs = VertexInfo(a-1,innerv,0)
            opstring(p) = 6*(b+1) + asubs - 1
        elseif ((r>W(0)/sum(W)) .and. (r<=(W(0)+W(1))/sum(W))) then
        ! the exit leg is chosen as 1 by Monte-Carlo
            if (W(1)==0) then
                write(*,*) 'wrong element'
            endif
            v2 = 4*p + 1
            asubs = VertexInfo(a-1,innerv,1)
            opstring(p) = 6*(b+1) + asubs - 1
        elseif ((r>(W(0)+W(1))/sum(W)) .and. (r<=(W(0)+W(1)+W(2))/sum(W))) then
        ! the exit leg is chosen as 2 by Monte-Carlo
            if (W(2)==0) then
                write(*,*) 'wrong element'
            endif
            v2 = 4*p + 2
            asubs = VertexInfo(a-1,innerv,2)
            opstring(p) = 6*(b+1) + asubs - 1
        elseif ((r>(W(0)+W(1)+W(2))/sum(W)) .and. (r<=1)) then
        ! the exit leg is chosen as 3 by Monte-Carlo
            if (W(3)==0) then
                write(*,*) 'wrong element'
            endif
            v2 = 4*p + 3
            asubs = VertexInfo(a-1,innerv,3)
            opstring(p) = 6*(b+1) + asubs - 1
        endif
        if (v1/=v2) loop_size=loop_size+1
        v1 = vertexlist(v2)
        call system_clock(int_time2)
        if (v1==v0) then
        ! the loop closed
            k = k+1
            exit
        elseif (loop_size>=100*mm) then
        ! the loop is too large before closing. Go back to the previous configuration
            opstring(:)=stringcopy(:)
            loop_size=0
            exit
        endif
    enddo
    total_loopsize=total_loopsize+loop_size
    if (k>=Nloop) then
        exit
    endif
enddo
 do i=0,nn-1
    ! Adjust the spin state according to new operator list after flipping loops
    ip = 0
        if (frstspinop(i) /= -1) then
            p = frstspinop(i)/4
            innerv = int(mod(frstspinop(i),4))
            a = mod(opstring(p),6) + 1
            b = opstring(p)/6 - 1
            if ((a /= 5) .and. (a /= 6)) then
                if (mod(a,2)==1) then
                    spin(i) = (2-a)*((-1)**innerv)
                else
                    spin(i) = (2-(mod(a,4)+1)) * ((-1)**innerv)
                endif
            elseif (a==5) then
                spin(i) = 1
            elseif (a==6) then
                spin(i) = -1
            endif
        else
            if (ran()<0.5) then
                spin(i)=-spin(i)
            endif
        endif
 enddo

end subroutine LoopUpdate
!------------------------!


!--------------------------!
subroutine MeasureWinding()        
!--------------------------!
! Measuring winding number !
use configuration;implicit none
integer :: p,op,b,a
W_Loop = 0
do p=0,mm-1
    op=opstring(p)
    b=op/6-1
    a=mod(op,6)+1
    if (b==nn-1) then
        if (a==2) then
            W_Loop=W_Loop-1
        elseif (a==4) then
            W_Loop=W_Loop+1
        endif
    endif
enddo

end subroutine MeasureWinding
!----------------------------!


!------------------------!
subroutine adjustcutoff(j)
!------------------------!
use configuration;implicit none
integer :: mnew,j
integer,allocatable :: stringcopy(:)

mnew = nh+nh/3

if (mnew>mm) then
    allocate(stringcopy(0:mm-1))
    stringcopy(:)=opstring(:)
    deallocate(opstring)
    allocate(opstring(0:mnew-1))
    opstring(0:mm-1)=stringcopy(:)
    opstring(mm:mnew-1)=0
    deallocate(stringcopy)
    mm=mnew
endif
deallocate(vertexlist)
allocate(vertexlist(0:4*mm-1))
vertexlist(:)=0

if (total_loopsize>1000) then
    Nloop=2000 * mm * Nloop / total_loopsize
    total_loopsize=0
endif

end subroutine adjustcutoff
!--------------------------!


!----------------------!
function IfMinus(a,inv,j)
!----------------------!
! Some invalid vertices (e.g. vertices that violate the spin conservation) !
! for the model are labeled by -1. This function is to check if for the    !
! given vertex, entrance leg and exit leg, the resulting a_new is -1.      !
use configuration;implicit none
integer :: a,inv,j
integer :: IfMinus

if (VertexInfo(a-1,inv,j)/=-1) then
    IfMinus = 1
else
    IfMinus = 0
endif

return
end function IfMinus
!-----------------!

!----------------------!
function IfOff(a,inv,j)
!----------------------!
! Check if for the given vertex, entrance leg and exit leg, the resulting !
! a_new is off-diagonal.                                                  !
use configuration;implicit none
integer :: a,inv,j
integer :: IfOff

if ((VertexInfo(a-1,inv,j)==2) .or. (VertexInfo(a-1,inv,j)==4)) then
    IfOff = 1
else
    IfOff = 0
endif

return
end function IfOff
!-----------------!




!-------------------------!
subroutine deallocateall()
!-------------------------!
use configuration;implicit none

deallocate(bsites)
deallocate(spin)
deallocate(opstring)
deallocate(frstspinop)
deallocate(lastspinop)
deallocate(vertexlist)

end subroutine deallocateall
!---------------------------!


 real function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
 implicit none

   !ran64=ran64*mul64+add64
   !ran=0.5d0+dmu64*dble(ran64)
 call random_seed()
 call random_number(ran)

 end function ran
!----------------!