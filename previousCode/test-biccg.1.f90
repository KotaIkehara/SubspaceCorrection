!implicit real*8(A-H,O-Z)  
integer,parameter :: maxnd=40000,maxns=1000000
integer,external :: blockiccg
character key*11
character mytype*3
character ver*6
character solver*128
character tag*32
character cdate*14
character text*128
integer ::  pr,nd,ns,ndb,ndx
double precision :: x0(maxnd),b(maxnd),ad(maxnd),au(maxns),abrs(1),res(maxnd)
integer :: ulnd(maxns),ulnt(maxns),nstp(1)
!integer :: nrowptr(maxnd+1)
integer :: i,rtc,lop1,lu0sw,mstp,setcoln,urptorg(maxnd+1)
double precision :: ganma,eps,maxb,rnorm,bnorm

      open(10,file='ieej.upper.qmf')
!      open(10,file='poison-crs.400.upper.qmf')
      read(10,*) key,mytype,ver
      read(10,*) solver
      read(10,*) tag
      read(10,*) cdate
      read(10,*) text
      read(10,*) pr
      read(10,*) nd,ns,ndb,ndx
      read(10,*)(ad(i),i=1,nd)
      read(10,*)(ulnd(i),i=1,nd)
      read(10,*)(au(i),i=1,ns)
      read(10,*)(ulnt(i),i=1,ns)
      if(ndb.gt.0)read(10,*)(b(i),i=1,nd)
      if(ndx.gt.0)read(10,*)(x0(i),i=1,nd)
      close(10)

      if (ndx.eq.0) then

        do i=1,nd
          x0(i)=0d0
		enddo

	  endif


      urptorg(1)=1
      do i=1,nd
        urptorg(i+1)=urptorg(i)+ulnd(i)
      enddo


      write(6,*) nd,ns,ndb,ndx,text

!      nrowptr(1)=1 
!      do i=1,nd
!        nrowptr(i+1)=nrowptr(i)+lnd(i)
!      enddo

	   
!    !$OMP PARALLEL
mstp=500
eps=1D-7
lop1=0
lu0sw=0
ganma=0.03d0
setcoln=0
rtc=blockiccg(x0,abrs,nstp,ad,au,b,ulnt,ulnd,nd,ns,mstp,eps,lop1,lu0sw,ganma,setcoln)
!call cg(ad,alu,lnt,lnd,nrowptr,nd,ns,x0,b,r)


!    !$OMP END PARALLEL

write(6,*) 'nstp(1)=',nstp(1)
write(6,*) 'error code=',rtc

do i=1,nd
  res(i)=ad(i)*x0(i) 
enddo

do i=1,nd
 do j=urptorg(i),urptorg(i+1)-1
   jj=ulnt(j)
   res(i)=res(i)+au(j)*x0(jj)
   res(jj)=res(jj)+au(j)*x0(i)

 enddo

enddo
maxb=0
do i=1,nd
res(i)=b(i)-res(i)
maxb=max(maxb,abs(res(i)))
enddo

bnorm=0d0
rnorm=0d0
do i=1,nd
 bnorm=bnorm+b(i)*b(i)
 rnorm=rnorm+res(i)*res(i)
enddo

write(6,*) 'Mugen norm=',maxb,'2norm=',dsqrt(rnorm/bnorm)

      END