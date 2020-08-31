subroutine biccg(ad,alu,lnt,nrowptr,au,ulnt,unrowptr,nd,ns,x,b,r,mstp,nstp,ganma,eps,ieer)
implicit none
integer :: i,nd,ns,ieer,lnt(*),nrowptr(*),ulnt(*),unrowptr(*)
integer :: interd,nstp(*),mstp,inum,ite
integer :: myid,idiv,omp_get_thread_num,omp_get_num_threads
integer :: numprocs,inps(0:100),inpe(0:100) 
double precision :: ad(*),au(*),alu(*)
double precision :: b(*),x(*),r(*),ganma,eps
double precision :: cgropp,cgrop,beta,alpha,bnorm,rnorm
double precision, dimension (:), allocatable :: pn,q
double precision, dimension (:), allocatable :: adic,bauic
double precision, dimension (:), allocatable :: z
integer, dimension(:), allocatable :: bnrowptr,blnt
external omp_get_thread_num,omp_get_num_threads

allocate(pn(nd),q(nd),adic(nd),bauic(ns),z(nd),blnt(ns),bnrowptr(nd+1),stat=ieer)
if (ieer.ne.0) then
ieer=900
return
endif


!$OMP PARALLEL PRIVATE(myid,idiv,inum,interd)
!!$OMP SINGLE
numprocs=omp_get_num_threads()
!!$OMP END SINGLE

 write(6,*) "numprocs=",numprocs

myid=omp_get_thread_num()

!$OMP BARRIER


!write(6,*) myid,numprocs

      interd=nd/numprocs+1
      interd=nd/numprocs
      do inum=0,numprocs-1
       inps(inum)=interd*inum+1
       inpe(inum)=interd*(inum+1)
      enddo
       inpe(numprocs-1)=nd

!write(6,*) "myid=",myid,inps(myid),inpe(myid),nd

call mkbu(ad,au,ulnt,unrowptr,adic,bauic,bnrowptr,blnt,nd,inps,inpe,myid,ganma)

call bic(adic,bnrowptr,blnt,bauic,inps,inpe,myid,ieer)

if (ieer.ne.0) then
mstp=0 
endif


!$OMP SINGLE
bnorm=0d0
!$OMP END SINGLE

!Calcurate right-hand vector norm
call pip(b,b,nd,bnorm) 

if (bnorm.lt.1d-15) then
  ieer=102
  mstp=0 
endif

! Set residula vector
call pmvp(ad,alu,lnt,nrowptr,nd,x,r)


!$OMP DO
do i=1,nd
r(i)=b(i)-r(i)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Enter iteration part
do ite=1,mstp
call fbsub(bnrowptr,blnt,bauic,nd,adic,myid,inps,inpe,z,r)


!$OMP SINGLE
cgropp=cgrop
cgrop=0d0
!$OMP END SINGLE

! !$OMP DO
!do i=1,nd
!  z(i)=r(i)
!enddo
!!$OMP BARRIER


call pip(r,z,nd,cgrop) 

!write(6,*) 'cgrop=',cgrop
if (ite==1) then

!$OMP DO
do i=1,nd
pn(i)=z(i)
enddo

else

!$OMP SINGLE
beta=cgrop/cgropp
!$OMP END SINGLE

!$OMP DO
do i=1,nd
  pn(i)=z(i)+beta*pn(i)
enddo

endif
call pmvp(ad,alu,lnt,nrowptr,nd,pn,q)

!$OMP SINGLE
alpha=0d0
!$OMP END SINGLE

call pip(pn,q,nd,alpha)


!$OMP SINGLE
alpha=cgrop/alpha
!$OMP END SINGLE

!write(6,*) 'alpha=',alpha,cgrop

!$OMP DO
do i=1,nd
  x(i)=x(i)+alpha*pn(i)
  r(i)=r(i)-alpha*q(i)
enddo


!$OMP SINGLE
rnorm=0d0
!$OMP END SINGLE
call pip(r,r,nd,rnorm)

!write(6,*) ite, dsqrt(rnorm/bnorm) 

if (dsqrt(rnorm/bnorm) < eps) then
  nstp(1)=ite
exit
endif
enddo

!$OMP END PARALLEL

if (ieer.eq.0.and.nstp(1).eq.0) then
nstp(1)=mstp
ieer=103
endif


deallocate(pn,q,adic,bauic,z,blnt,bnrowptr)


return
end





subroutine fbsub(bnrowptr,blnt,bauic,nd,adic,myid,inps,inpe,z,r)
integer :: i,j,inps(0:*),inpe(0:*),myid
integer :: blnt(*),bnrowptr(*)
double precision :: adic(*),bauic(*),z(*),r(*)

!FORWARD SUBSTITUTION
!$OMP DO
do i=1,nd
  z(i)=r(i)
enddo

do i=inps(myid),inpe(myid)-1
  do j=bnrowptr(i),bnrowptr(i+1)-1
    jj=blnt(j)
    z(jj)=z(jj)-z(i)*bauic(j)/adic(i)
  enddo
enddo

!BACKWARD SUBSTITUTION
z(inpe(myid))=z(inpe(myid))/adic(inpe(myid))

do i=inpe(myid)-1,inps(myid),-1
  do j=bnrowptr(i),bnrowptr(i+1)-1
    jj=blnt(j)
    z(i)=z(i)-bauic(j)*z(jj)
  enddo
    z(i)=z(i)/adic(i)
enddo



!$OMP BARRIER 


return
end


subroutine mkbu(ad,au,ulnt,unrowptr,adic,bauic,bnrowptr,blnt,nd,inps,inpe,myid,ganma)
implicit none
integer :: i,j,kk,jstart,inps(0:*),inpe(0:*),myid,nd
integer :: ulnt(*),unrowptr(*),blnt(*),bnrowptr(*)
double precision :: ad(*),au(*),adic(*),bauic(*),ganma

!write(6,*) "myid=",myid,"inps=",inps(myid),"inpe=",inpe(myid),ganma

do i=inps(myid),inpe(myid)
  adic(i)=ad(i)*(ganma+1d0)
enddo

kk=0
jstart=unrowptr(inps(myid))
bnrowptr(inps(myid))=jstart
do i=inps(myid),inpe(myid)-1
  do j=unrowptr(i),unrowptr(i+1)-1
    if (ulnt(j).le.inpe(myid)) then
      blnt(kk+jstart)=ulnt(j)
      bauic(kk+jstart)=au(j)
      kk=kk+1
    endif
enddo
 bnrowptr(i+1)=kk+jstart
enddo

!$OMP BARRIER
!write(6,*) myid,bnrowptr(inps(myid)), bnrowptr(inpe(myid))-1

return
end

subroutine bic(adic,bnrowptr,blnt,bauic,inps,inpe,myid,ieer)
implicit none
integer :: i,j,jp,jjp,jj,ji,ielocal,ieer
integer :: myid,inpe(0:*),inps(0:*)
integer :: bnrowptr(*),blnt(*)
double precision :: adic(*),bauic(*) 

ielocal=0

do i=inps(myid),inpe(myid)-1
  if (ielocal.eq.0) then
  do j=bnrowptr(i),bnrowptr(i+1)-1
    jj=blnt(j)
    adic(jj)=adic(jj)-bauic(j)*bauic(j)/adic(i)

      if (adic(jj).lt.1d-8) then
        ielocal=101
      endif

    do jp=bnrowptr(i),bnrowptr(i+1)-1
      jjp=blnt(jp)
      if (jjp > jj) then
        do ji=bnrowptr(jj),bnrowptr(jj+1)-1
          if (blnt(ji).eq.jjp) then
            bauic(ji)=bauic(ji)-bauic(j)*bauic(jp)/adic(i)
            exit
          endif
        enddo
      endif
    enddo
  enddo
  endif
enddo
!$OMP BARRIER
 
if (ielocal.ne.0.and.ieer.eq.0) ieer=ielocal

return
end