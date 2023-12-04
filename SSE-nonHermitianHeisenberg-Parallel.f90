!--------------------!
 module configuration
!--------------------!
 save

 integer :: nn       ! number of sites;
 integer :: nb       ! number of bonds; nb=2*nn
 integer :: if_PBC

 integer,parameter :: msteps=20000
 integer,parameter :: nbins=3
 integer,parameter :: isteps=80000

 integer :: nh      
 integer :: mm
 integer :: Nloop
 integer :: W_Loop=0
 integer :: W_Loop2=0
 integer :: total_loopsize=0
 integer :: N_WSector(0:nbins-1,-40:40)
 integer :: N_spinsector(0:nbins-1,-10:10)

 integer, allocatable :: bsites(:,:)

 integer, allocatable :: frstspinop(:) ! first operation on each site in linked vertex list
 integer, allocatable :: lastspinop(:)
 integer, allocatable :: vertexlist(:)

 integer, allocatable :: spin(:)
 integer, allocatable :: opstring(:)

 real(8),dimension(0:5,0:3,0:3) :: VertexInfo
 real(8),dimension(0:5) :: MatrixElement

 real(8) :: Jz
 real(8) :: dJ
 real(8) :: delta
 real(8) :: beta
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

 real(8) :: dJ_out,delta_out,E
 real(8) :: S,EL(0:nbins-1),EAndVar(2)
!  real(8) :: EList(0:2,-7:7,0:nbins-1)
!  real(8) :: EList_FF(0:2,-7:7,0:nbins-1)
!  real(8) :: WindingPeakList(0:2,-7:7)
 character*8 :: str

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,p,ierr)
call MPI_Barrier(MPI_COMM_WORLD,ierr)

write(str,'(I4)') my_rank

open(unit=10,file='Results_PBC_varJz'//trim(adjustl(str))//'.txt',status='replace')


dJ_out=0.1d0*(my_rank/15)+0.3
delta_out=0.05d0*mod(my_rank,15)
! Jz_out=0
! delta_out=0.7*dble((-1)**(my_rank))

call main(dJ_out,delta_out,EL)

! do j=0,2
!     dJ=dble(j)/dble(10)
!     do i=0,7
!         delta=dble(i)/dble(10)
!         call main(EL,wp)
!         EList(j,i,:)=EL(:)
!         WindingPeakList(j,i)=wp
        write(10,*) 'Jz=',Jz,'  dJ=',dJ,'  delta=',delta,'  -E/N=',EAndVar(1),'  variation:',EAndVar(2)
        write(10,*) 'wp=',wp
        write(10,*) N_WSector
!     enddo
! enddo

close(10)
call MPI_Barrier(MPI_COMM_WORLD,ierr)

call MPI_FINALIZE(ierr)
 end program SSERun
!=====================!


!----------------!
subroutine main(dJ_out,delta_out,EAndVar)
!----------------!
use configuration;implicit none
integer :: i,j,k,wp,wpL(1)
real(8) :: dJ_out,delta_out
real(8) :: E,Eavg,EAndVar(2),EList(nbins)

N_WSector(:,:)=0
N_spinsector(:,:)=0
nn=64
Jz=0
delta=delta_out
dJ=dJ_out
beta=100
! delta=0.5
if_PBC=1

MatrixElement=(/ C+Jz*0.25d0, 0.5d0*(1-delta), C+Jz*0.25d0, 0.5d0*(1+delta), C-Jz*0.25d0, C-Jz*0.25d0/)
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
        if (abs(sum(spin))<=10) then
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

write(*,*) 'Jz=',Jz,'  dJ=',dJ,'  delta=',delta,'  -E/N=',EAndVar(1),'  variation:',EAndVar(2)
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
use configuration;implicit none
integer,dimension(0:3) :: vertexspin,vertexspin_aft
integer :: a,xin,xout

do a=0,5
    if ((a==0) .or. (a==3) .or. (a==4)) then
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

!   do i=0,nn-1
!     spin(i)=(-1)**(i+1)
!  enddo

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
use configuration;implicit none
integer,dimension(0:nn-1) :: spin0_1
integer :: p,i_op,op,b,a,i,ip,j_A
real(8) :: rd,Pa,Pd,aprob,dprob

spin0_1(:)=spin(:)
do p=0,mm-1
    op = opstring(p)
    if (op==0) then
        b = int(ran()*nb)
        rd = ran()
        Pa = rd*(mm-nh)
        if (Pa < nb*beta*(C-spin(bsites(0,b))*spin(bsites(1,b))*Jz*0.25d0)) then
            if (spin(int(bsites(0,b)))==1) then
                if (spin(int(bsites(1,b)))==-1) then
                    a = 1
                elseif (spin(int(bsites(1,b)))==1) then
                    a = 5                        
                endif
            elseif (spin(int(bsites(0,b)))==-1) then
                if (spin(int(bsites(1,b)))==-1) then
                    a = 6
                elseif (spin(int(bsites(1,b)))==1) then
                    a = 3
                endif
            endif
            opstring(p)=6*(b+1) + a - 1
            nh = nh+1
            if ((spin(int(bsites(0,b)))*spin(int(bsites(1,b)))>0) .and. ((a==1) .or. (a==3))) then
                write(*,*) 'invalid operator'
            elseif ((spin(int(bsites(0,b)))*spin(int(bsites(1,b)))<0) .and. ((a==5) .or. (a==6))) then
                write(*,*) 'invalid operator'
            endif
        endif
    elseif (mod(op,2)==0 .or. (int(mod(op,6)+1)==6)) then
        b = int(op/6 -1)
        a = int(mod(op,6)+1)
        rd = ran()
        Pd = rd/dble(mm-nh+1)
        aprob = nb*beta*(C-spin(bsites(0,b))*spin(bsites(1,b))*Jz*0.25d0)
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
        write(*,*) 'time period is broken'
        exit
    endif
enddo
end subroutine diagonalupdate
!----------------------------!


!------------------------!
subroutine LinkedVertex()
!------------------------!
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
    stringcopy(:)=opstring(:)
    call system_clock(int_time1)
    loop_size=0
    do
        if (opstring(v0/4)==0) then
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
            if (W(0)==0) then
                write(*,*) 'wrong element'
            endif
            v2 = 4*p + 0
            asubs = VertexInfo(a-1,innerv,0)
            opstring(p) = 6*(b+1) + asubs - 1
        elseif ((r>W(0)/sum(W)) .and. (r<=(W(0)+W(1))/sum(W))) then
            if (W(1)==0) then
                write(*,*) 'wrong element'
            endif
            v2 = 4*p + 1
            asubs = VertexInfo(a-1,innerv,1)
            opstring(p) = 6*(b+1) + asubs - 1
        elseif ((r>(W(0)+W(1))/sum(W)) .and. (r<=(W(0)+W(1)+W(2))/sum(W))) then
            if (W(2)==0) then
                write(*,*) 'wrong element'
            endif
            v2 = 4*p + 2
            asubs = VertexInfo(a-1,innerv,2)
            opstring(p) = 6*(b+1) + asubs - 1
        elseif ((r>(W(0)+W(1)+W(2))/sum(W)) .and. (r<=1)) then
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
            k = k+1
            exit
        elseif (loop_size>=100*mm) then
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


!------------------------------!
subroutine LoopUpd_Determined()
!------------------------------!
use configuration;implicit none

integer :: i,v0,v1,v2,looplen,j
 integer :: totallooplen=0
 integer :: W_singleLoop=0
 integer :: W_loopCopy,W_loopCopy2
 integer :: p,b,a,innerv,asubs,varW
 integer, allocatable :: stringcopy(:)
 integer, allocatable :: v1list(:)
 integer, allocatable :: v2list(:)
 integer, allocatable :: iv1list(:)
 integer, allocatable :: iv2list(:)
 real :: r,ran
 real :: Ploop,W1,W2
 
 W_loopCopy=W_Loop
 W_loopCopy2=W_loopCopy
 allocate(stringcopy(0:mm-1))
 W_Loop2=0
 
 do v0=0,4*mm-1,2       
    if (vertexlist(v0)<0) cycle                  
    v1=v0
    stringcopy(:)=opstring(:)
    totallooplen=totallooplen+looplen*2
    Ploop = 1.d0/dble(4*mm)
    looplen=0                                                               
    do                                    
      innerv=mod(v1,4)
      p=v1/4
      b=stringcopy(p)/4
      a=mod(stringcopy(p),4)+1
      if (innerv/2==0) then
         asubs=ieor(a-1,1)+1
      else
         asubs=5-a
      endif
      W1=MatrixElement(asubs)-0.5*dJ*dble(mod(asubs-1,2))*dble((-1)**(bsites(1,b)-1))
      W2=MatrixElement(a)-0.5*dJ*dble(mod(a-1,2))*dble((-1)**(bsites(1,b)-1))
      Ploop=Ploop*(W1/(W1+W2))
      stringcopy(p)=4*b+asubs-1
      if (b==nn) then
         ! write(*,*) 'PBC boundary crossed'
         if (asubs==2) then
            W_loopCopy = W_loopCopy-1
         elseif (asubs==4) then
            W_loopCopy = W_loopCopy+1
         elseif (a==2) then
            W_loopCopy = W_loopCopy+1
         elseif (a==4) then
            W_loopCopy = W_loopCopy-1
         endif
      endif
      v2=ieor(v1,1)
      if (v2>v1) then
         W_Loop2=W_Loop2+1
      else
         W_Loop2=W_Loop2-1
      endif

      if (looplen==0) then
         allocate(v1list(1))
         allocate(v2list(1))
         v1list(1)=v1
         v2list(1)=v2
         looplen=looplen+1
      else
         allocate(iv1list(looplen))
         allocate(iv2list(looplen))
         iv1list(:)=v1list(:)
         iv2list(:)=v2list(:)
         deallocate(v1list)
         deallocate(v2list)
         allocate(v1list(looplen+1))
         allocate(v2list(looplen+1))
         v1list(1:looplen)=iv1list(:)
         v1list(looplen+1)=v1
         v2list(1:looplen)=iv2list(:)
         v2list(looplen+1)=v2
         deallocate(iv1list)
         deallocate(iv2list)
         looplen=looplen+1
      endif

      v1=vertexlist(v2)

      if (v1==v0) exit
    enddo
    ! TotalLoop=TotalLoop+1
    ! if (W_loopCopy>W_loopCopy2) then
    !   UpGlobalLoop=UpGlobalLoop+1
    !   GlobalLoop=GlobalLoop+1
    ! elseif (W_loopCopy<W_loopCopy2) then
    !   DownGlobalLoop=DownGlobalLoop+1
    !   GlobalLoop=GlobalLoop+1
    ! endif
    
    Ploop = Ploop*4*dble(mm)*dble(2**(looplen-2))

      ! open(33,file='Ploop.txt',status='replace')
      ! write(33,*) Ploop,log(Ploop)
      ! close(33)
      ! Ploop = Ploop*4*mm*(((0.05+0.5*(1-delta))/(0.5*(1-delta))*((0.05+0.5*(1+delta))/(0.5*(1+delta))))**(looplen/2-2))

    ! if (N_totalLoop==1) then
    !   GLoopLen=looplen
    !   MLoopLen=looplen
    !   minPloop=Ploop
    !   maxPloop=Ploop
    ! else
    !   if (Ploop<minPloop) then
    !      minPloop=Ploop
    !      MLoopLen=looplen
    !   elseif (Ploop>maxPloop) then
    !      maxPloop=Ploop
    !      GLoopLen=looplen
    !   endif
    ! endif
    ! N_totalLoop=N_totalLoop+1

   !  if (b==nb) then
   !    write(*,*) Ploop
   !  endif
    r=ran()
    write(18,*) r
    ! if (W_loopCopy2/=W_loopCopy) then
    !   varW=(int(sign(1,W_loopCopy-W_loopCopy2))+1)/2
    !   Rate_WSector(varW,W_loopCopy2)=Rate_WSector(varW,W_loopCopy2)+Ploop
    !   NGlobalLoop_WSector(varW,W_loopCopy2)=NGlobalLoop_WSector(varW,W_loopCopy2)+1
    !   ! write(44,*) 'W_before=',W_loopCopy2,' W_after=',W_loopCopy,' r=',r,' Ploop=',Ploop,' NGlobalUp=',NGlobalLoop_WSector(1,W_loopCopy2),' NGlobalDown=',NGlobalLoop_WSector(0,W_loopCopy2)
    ! endif
    ! NTotalLoop_WSector(W_loopCopy2)=NTotalLoop_WSector(W_loopCopy2)+1
    if (r<Ploop) then
      do i=1,looplen
         v1=v1list(i)
         v2=v2list(i)
         vertexlist(v1)=-2
         vertexlist(v2)=-2
         opstring(:)=stringcopy(:)
      enddo
    else
      do i=1,looplen
         v1=v1list(i)
         v2=v2list(i)
         vertexlist(v1)=-1
         vertexlist(v2)=-1
      enddo
    endif
    W_loopCopy2=W_loopCopy
    deallocate(v1list)
    deallocate(v2list)
 enddo

 do i=1,nn
    if (frstspinop(i)/=-1) then
       if (vertexlist(frstspinop(i))==-2) spin(i)=-spin(i)  !表示site i in |α(0)>属于之前经过了算符替换，且穿过了time period的loop,所以i的自旋也进行翻转!
    else
       if (ran()<0.5) then
         spin(i)=-spin(i)
       endif
    endif
 enddo
!  write(18,*) '   '
!  W_Loop2=W_Loop2/nn
!  write(22,*) W_Loop2
!  write(22,*) '   '
!  write(44,*) '----------------------------'
 deallocate(stringcopy)

end subroutine LoopUpd_Determined
!--------------------------------!


!--------------------------!
subroutine MeasureWinding()
!--------------------------!
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