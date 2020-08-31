integer*4 function blockiccg(x,abrs,nstp,ad,au,b,ulnt,ulnd,nd,ns,mstp,eps,lop1,lu0sw,ganma,setcoln)

implicit none
integer :: i,nd,ns,ieer,ulnt(*),ulnd(*)
integer :: nstp(*),mstp,lop1,lu0sw,setcoln

double precision :: ad(*),au(*),abrs(*),ganma
double precision :: b(*),x(*),eps
double precision, dimension (:), allocatable :: r

integer, dimension (:), allocatable :: llnd,lnrowptr,llnt,unrowptr
double precision, dimension (:), allocatable :: al

integer, dimension (:), allocatable :: lnd,nrowptr,lnt
double precision, dimension (:), allocatable :: alu


setcoln=0
write(6,*) 'ENTER BLOCK-ICCG'



allocate(unrowptr(nd+1),stat=ieer)
if (ieer.ne.0) then
blockiccg=900
return
endif

unrowptr(1)=1 
do i=1,nd
  unrowptr(i+1)=unrowptr(i)+ulnd(i)
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Generate lower-triangular part of coefficient matrix 

allocate(llnd(nd),lnrowptr(nd+1),llnt(ns),al(ns),stat=ieer)
if (ieer.ne.0) then
blockiccg=900
return
endif

call genlow(au,ulnt,unrowptr,nd,al,llnt,llnd,lnrowptr,ieer)

if (ieer.ne.0) then
blockiccg=ieer
return
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check symmetricy 
!
!do i=1,nd
!  do j=lnrowptr(i),lnrowptr(i+1)-1
!    jj=llnt(j)
!    kk=0
!	do ji=unrowptr(jj),unrowptr(jj+1)-1

!      if (ulnt(ji).eq.i) then
!        if (al(j).ne.au(ji)) then
!        write(6,*) "error not match",i,jj,al(j),au(ji)
!        endif
!        kk=1
!	  endif
!	enddo
!	   if (kk.eq.0) write(6,*) "error not found",i,jj
    
!  enddo
!enddo

!write(6,*) lnrowptr(nd+1)-1,ns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generate all part of coefficient matrix

allocate(lnd(nd),nrowptr(nd+1),lnt(ns*2),alu(ns*2),stat=ieer)
if (ieer.ne.0) then
blockiccg=900
return
endif

call genall(au,ulnt,unrowptr,nd,al,llnt,lnrowptr,alu,lnt,nrowptr,lnd,ieer)

if (ieer.ne.0) then
blockiccg=ieer
return
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate array for low part
deallocate(llnt)
deallocate(lnrowptr)
deallocate(llnd)
deallocate(al)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


write(6,*) "Setup is finished" 

allocate(r(nd),stat=ieer)
if (ieer.ne.0) then
blockiccg=900
return
endif

call biccg(ad,alu,lnt,nrowptr,au,ulnt,unrowptr,nd,ns,x,b,r,mstp,nstp,ganma,eps,ieer)

deallocate(lnd,nrowptr,lnt,alu)
deallocate(r)



blockiccg=ieer

return
end function blockiccg 