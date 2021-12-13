program voffunc
use vof2d
use dflib
implicit none
integer L1,M1,L2,M2,L3,M3,i,j
real(8) xmin,xmax,ymin,ymax,fs(ni,nj),x(ni),y(nj),dx,dy,cfl,dt,t,tend,ui,vi
real(8) snormx(ni,nj),snormy(ni,nj),u(ni,nj),v(ni,nj),p(ni,nj),ax(ni,nj),ay(ni,nj),fb(ni,nj)
real(8) nx(ni,nj),ny(ni,nj),c(ni,nj),EL1,EL2,ct,cnw,cne,csw,cse,vol(ni,nj) !,v0
integer num,itimestep,prt, mark(ni,nj)
real(8) fn(ni,nj), init_f,tmp,tmp1
logical,external:: domain
logical tr
character*200 filename1,filename2
integer(2)::control,clearcontrol,newcontrol

call getcontrolfpqq(control)
clearcontrol=control.and.(.not.fpcw$mcw_pc)
newcontrol=clearcontrol.or.fpcw$64
call setcontrolfpqq(newcontrol)



xmin=0d0;xmax=1d0;ymin=0d0;ymax=1d0;dt=0.025;t=0d0;tend=1d0
num=0;prt=1;cfl=0.25d0;
ax=1d0;ay=1d0;fb=1d0;fs=0d0
EL1=0d0;EL2=0d0
!open( 10,file="rezult.dat",status='new')
call setup()
call grid(xmin,xmax,ymin,ymax,x,y,dx,dy,L1,M1)
L2=L1-1;M2=M1-1; L3=L2-1;M3=M2-1
call fss(x,y,dx,dy,fs,L1,M1)
do i=1,L1
  do j=1,M1
    fn(i,j)=fs(i,j)
    init_f = init_f+fs(i,j) 
  end do
end do

call marking(fb,fs,mark,L1,M1)
call normal(.true.,L1,M1,x,y,fs,mark,snormx,snormy)
do i=2,L1-1
    do j=2,M1-1
	 if(mark(i,j)==S) then
       call planec(fs(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
	 end if
	 if(fs(i,j)<=0d0+1d-5) then	    
		
		vol(i,j)=0d0
	  else if(fs(i,j)>=1d0-1d-5)then
	   
		vol(i,j)=dx*dy
	  else
	    if(mark(i,j)/=E.or.mark(i,j)/=F)then
		    
			call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse,vol(i,j))
			
			if(vol(i,j)<0d0) vol(i,j)=0d0
		  		    
		end if
		  			
	  end if
	end do
end do

do i=2,L1-1
 do j=2,M1-1
   vnull=vnull+vol(i,j)
 end do
end do



filename1=Trim("../voffunc3/vtk/")//trim('uvpf_null')//trim('.vtk')


filename2=Trim("../voffunc3/vtk/")//trim('surf_null')//trim('.vtk')

call vof_vtk(100,filename1,filename2,u,v,p,fs,mark,x,y,snormx,snormy,c,L1,M1,dx,dy)





do i=1,L1
  call uvp(4,x(i),y(1),ui,vi)
  u(i,1)=ui
  v(i,1)=vi
end do

do j=1,M1
  call uvp(4,x(1),y(j),ui,vi)
  u(1,j)=ui
  v(1,j)=vi
end do

do i=2,L1
 do j=2,M1
   call uvp(4,x(i),y(j),ui,vi)
   u(i,j)=ui
   v(i,j)=vi
 end do
end do
!!!!!!!!!!!!!!!test
call plc(-0.707d0,0.707d0,0.0d0,10.0d0,10.0d0,0.0d0,cnw,csw,cne,cse)
call value2(3.55d0,-0.707d0,0.707d0,0.0d0,10.0d0,10.0d0,0.0d0,cnw,csw,cne,cse,vol)
    !value2(c,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,v)
!print*,'vol=', vol
!call planec(fs1,nx,ny,xw,xe,yn,ys,c)
call planec(0.876d0,-0.707d0,0.707d0,0.0d0,10.0d0,10.0d0,0.0d0,ct)
!print*,ct

!!!!!!!!!!!!!!!!!!!!!!!!!!!
itimestep=0
num=1
do while(t<tend)
  call curant(u,v,L1,M1,dx,dy,dt,cfl,t)    
  !call advect_da4(u,v,fs,fb,mark,ax,ay,snormx,snormy,c,x,y,L1,M1,dx,dy,dt,t,itimestep)
  call advect_da(u,v,fs,fb,mark,ax,ay,snormx,snormy,c,x,y,L1,M1,dx,dy,dt,t,itimestep)
    
if(mod(num,prt)==0) then     
  call output(1,u,v,p,fs,mark,x,y,snormx,snormy,c,L1,M1,dx,dy,num)
end if

tmp=0d0
do i=1,L1
  do j=1,M1    
    tmp = tmp+fs(i,j) 
  end do
end do

write(*,'(1x,a,i5,1x,a,d10.5,1x,a,d40.33)') 'step=',itimestep, 't=', t,  'mass=',tmp/init_f 
  t=t+dt
  num=num+1
  itimestep=itimestep+1  
end do

tmp=0.0d0
tmp1=0.0d0
do j=1,L1
 do i=1,M1
   tmp = tmp + dabs(fs(i,j)-fn(i,j))
   tmp1 = tmp1 + fs(i,j)
 end do
 end do
tmp=tmp/init_f
tmp1 = abs((tmp1-init_f))/init_f
write(*,*) '# Error L1 norm=',tmp 
write(*,*) '# Mass error=',tmp1 


!call PRNT(fs,"vof_fraction",10,L1,M1)
end

subroutine grid(xmin,xmax,ymin,ymax,x,y,dx,dy,L1,M1)
use vof2d
real(8),intent(in)::xmin,xmax,ymin,ymax
real(8),intent(out)::x(ni),y(nj),dx,dy
integer,intent(out)::L1,M1
integer L2,M2,L3,M3,i,j
L1=101;M1=101                         !L1=13;M1=13
L2=L1-1;M2=M1-1;L3=L1-2;M3=M1-2
dx=(xmax-xmin)/float(L3)
dy=(ymax-ymin)/float(M3)

do i=1,L2
  x(i)=xmin+(i-1)*dx
end do

do j=1,M2
  y(j)=ymin+(j-1)*dy
end do

end subroutine grid

logical function domain2(x,y)
real(8),intent(in):: x,y
real(8) a,b,R
a=0.5d0; b=0.5d0;  R=0.1d0                     ! R=0.25d0
domain2=((x-a)*(x-a)+(y-b)*(y-b)<=R*R)
return
end function domain2

logical function domain(x,y)
real(8),intent(in):: x,y
real(8) a,b,R,c,xe,xw,ys,yn
logical circle,rectangle
a=0.5d0; b=0.75d0; R=0.15d0; c=1d0/3d0*R    !R=0.15d0
!circle=((x-a)*(x-a)+(y-b)*(y-b)<=R*R)
!xe=a+c/2d0; xw=a-c/2d0;ys=b-R;yn=b+(R-c)
!rectangle=(xw<x<xe).and.(ys<y<yn)
!domain=circle.xor.rectangle
domain=sqrt((x-0.5d0)**2+(y-0.75d0)**2)<=0.15.and.(abs(x-0.5d0)>0.03d0.or.y>0.85d0)

return
end function domain

subroutine fss(x,y,dx,dy,fs,L1,M1)
use vof2d
real(8),intent(in):: x(ni),y(nj),dx,dy
integer,intent(in):: L1,M1
real(8),intent(out):: fs(ni,nj)
integer i,j,L2,M2,L3,M3
real(8) xne,yne,xse,yse,xnw,ynw,xsw,ysw,volume
real(8) dxl,dyl,vof,xd(nd),yd(nd)
logical mne,mse,mnw,msw
logical,external:: domain
L2=L1-1;M2=M1-1;L3=L1-2;M3=M3-2
fs=0d0
do i=2,L2
  do j=2,M2
    vof=0d0
    xne=x(i);   yne=y(j)
	xse=x(i);   yse=y(j-1)
	xnw=x(i-1); ynw=y(j)
	xsw=x(i-1); ysw=y(j-1)

	mne=domain(xne,yne)
	mse=domain(xse,yse)
	mnw=domain(xnw,ynw)
	msw=domain(xsw,ysw)

    volume=(x(i)-x(i-1))*(y(j)-y(j-1))

	if(mne.or.mse.or.mnw.or.msw) then
	  dxl=dx/nd
	  dyl=dy/nd
      xd(1)=x(i-1)+dxl/2d0
	  yd(1)=y(j-1)+dyl/2d0
      do k1=2,nd
	    xd(k1)=xd(1)+(k1-1)*dxl
	  end do
	  do k2=2,nd		  
	    yd(k2)=yd(1)+(k2-1)*dyl
	  end do

      
      

	  do k1=1,nd
	    do k2=1,nd		  
		  if(domain(xd(k1),yd(k2))) then
		    vof=vof+1d0/float(nd*nd)
		  end if
	    end do
	  end do
	end if
    fs(i,j)=vof
  end do
end do

end subroutine fss


subroutine fss2(x,y,dx,dy,fs,NX,NY)
use vof2d
integer,intent(in):: NX,NY
real(8),intent(in)::x(ni),y(nj),dx,dy
real(8),intent(out):: fs(ni,nj) 
real(8) tmp
integer ii
ii=10
do j=1,NY*ii
      do i=1,NX*ii
         tmp=(float(i)*dx/float(ii)-0.51d0)**2&
             +(float(j)*dy/float(ii)-0.76d0)**2
         if((sqrt(tmp).lt.0.17d0)&
            .and.((abs(float(i)*dx/float(ii)-0.51d0).gt.0.03d0)&
             .or.(float(j)*dx/float(ii).ge.0.85d0))) then  
            fs(i/ii+1,j/ii+1)=fs(i/ii+1,j/ii+1)+1.0d0/float(ii*ii)
         end if
      end do
      end do 

end subroutine fss2





SUBROUTINE PRNT(FI,TITLE,numfile,NNI,NNJ)
use vof2d
real(8),intent(in)::FI(ni,nj)    
CHARACTER*6,intent(in):: TITLE
integer,intent(in)::NNI,NNJ
      WRITE(numfile,20) TITLE 
      IS=-11
  100 IS=IS+12
      IE=IS+11
      IE=MIN(NNI,IE)
      WRITE(numfile,21) (I,I=IS,IE) 
      WRITE(numfile,22) 
      DO J=NNJ,1,-1
        WRITE(numfile,23) J,(FI(I,J),I=IS,IE) 
      END DO
      IF(IE.LT.NNI) GO TO 100
   20 FORMAT(2X,26('*-'),7X,A6,7X,26('-*')) 
   21 FORMAT(3X,'I = ',I3,11I10)
   22 FORMAT(2X,'J')
   23 FORMAT(1X,I3,1P12E10.2) 
      RETURN
      END
	  
	  
	  
	  
subroutine uvp(flag,x,y,ui,vi)
use vof2d
implicit none
integer,intent(in)::flag
real(8),intent(in)::x,y
real(8),intent(out)::ui,vi
real(8) omega,x0,y0
!1-diagonal,2-single vortex,3-deformation field
if(flag==1)then
  ui=1d0
  vi=0d0
else if(flag==2)then
  ui=-2d0*sin(pi*x)*sin(pi*x)*sin(pi*y)*cos(pi*y)
  vi=2d0*sin(pi*x)*sin(pi*y)*sin(pi*y)*cos(pi*x)
else if(flag==3)then
  ui=-(-sin(4d0*pi*(x+0.5d0)))*sin(4d0*pi*(y+0.5d0))
  vi=sin(4d0*pi*(y+0.5d0))*cos(4d0*pi*(y+0.5d0))
else if(flag==4) then
  omega=1d0;x0=0.5d0;y0=0.75d0
  
  ui=-omega*(y-y0)
  vi=omega*(x-x0)
end if
end subroutine  uvp
  