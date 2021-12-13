subroutine advect2(u,v,fs,fb,mark,ax,ay,x,y,L1,M1,dx,dy,dt,t,itimestep)
!slic donor acceptor method x-y sweep: u(i,j)<0 fs<0 ????
use vof2d
real(8) fc,fss
!integer,parameter::E=0,F=1,S=2
 
integer,intent(in):: L1,M1,itimestep
real(8),intent(in):: t,dt,u(ni,nj),v(ni,nj),x(ni),y(nj),dx,dy,fb(ni,nj),ax(ni,nj),ay(ni,nj)
real(8),intent(inout)::fs(ni,nj)
integer,intent(inout)::mark(ni,nj)
real(8) flux_x(ni,nj),flux_y(ni,nj),dvx(ni,nj),dvy(ni,nj),VL(ni,nj),vx(ni,nj),vy(ni,nj)
real(8) fsx(ni,nj),fsy(ni,nj)
real(8) fd,fdm1,fdp1,v0,dv
logical us
open(2,file="result.dat")


!call advect_da2(u,v,fs,fb,mark,ax,ay,x,y,L1,M1,dx,dy,dt,t,itimestep)


fc=0.01
fss=0.4
!dx=x(2)-x(1)
!dy=y(2)-y(1)
v0=dx*dy
flux_x=0d0
flux_y=0d0
fsx=0d0
fsy=0d0

if(mod(itimestep,2)==0) then
   write(2,*) 'x-y sweep'
   write(2,*) 'solving fluxes for x-axis'
   !solving fluxes for x-axis
do i=2,L1-1
 do j=2,M1-1
   
  
   if(u(i,j)>0d0) then
     fd=fs(i,j)
	 fdm1=fs(i-1,j)
	 fdp1=fs(i+1,j)	
   end if

   if(u(i,j)<0d0) then
     fd=fs(i,j)
	 fdm1=fs(i+1,j)
	 fdp1=fs(i-1,j)	
   end if
   
   !if(u(i,j)==0d0) exit



   if (fdm1<=fc.and.fd>=fc.and.fdp1>=fc) then
       flux_x(i,j)=min(fd*dx,abs(u(i,j))*dt*fdp1)
        
   end if

    if (fdm1>=1d0-fc.and.fd>=fc.and.fdp1<=fc) then
        flux_x(i,j)=min(fd*dx,max(0d0,abs(u(i,j))*dt-(1-fd)*dx))
		
   end if

    if (fdm1>=fc.and.fd>=fc.and.fdp1>=fc) then
        flux_x(i,j)=abs(u(i,j))*dt*fd
		
   end if

    if (fdm1>=1-fc.and.fd>=fc.and.fdp1>=fc) then
        flux_x(i,j)=min(fd*dx,abs(u(i,j))*dt*fdp1+max((1-fdp1)*abs(u(i,j))*dt-(1d0-fd)*dx,0d0))
		
   end if

    if (fdm1>=fc.and.fd>=fc.and.fdp1>=1-fc) then
       flux_x(i,j)=min(fdp1*abs(u(i,j))*dt, fd*dx-fdm1*(dx-abs(u(i,j))*dt))
	   
   end if
    !if(fdm1==0d0.and.fd>0d0.and.fdp1==0d0) then
	 !  flux_x(i,j)=min(fd*dx,abs(u(i,j))*dt*fd)
	!end if
    us=(fdm1<fc.and.fd<fc.and.fdp1<fc).or.(fdm1>fc.and.fd>fc.and.fdp1>fc).or.&
	   (fdm1<fc.and.fd<fc).or.(fd<fc.and.fdp1<fc)
	if(us) then
	   flux_x(i,j)=abs(u(i,j))*dt*fd
	end if
	 
 end do
end do

do i=2,L1-1
 do j=2,M1-1
   fsx(i,j)=(fs(i,j)+1d0/dx*(sign2(u(i-1,j))*flux_x(i-1,j)-sign2(u(i,j))*flux_x(i,j)))/(1d0-dt/dx*(u(i,j)-u(i-1,j)))

 end do 
end do 

do i=2,L1-1
 do j=2,M1-1  
  !solving fluxes for y-axis
  write(2,*) 'solving fluxes for y-axis'
   if(v(i,j)>0d0) then
     fd=fsx(i,j)
	 fdm1=fsx(i,j-1)
	 fdp1=fsx(i,j+1)
    
   end if

   if(v(i,j)<0d0) then
     fd=fsx(i,j)
	 fdm1=fsx(i,j+1)
	 fdp1=fsx(i,j-1)
     
   end if

   !if(v(i,j)==0d0) exit

   if (fdm1<=fc.and.fd>=fc.and.fdp1>=fc) then
       flux_y(i,j)=min(fd*dy,abs(v(i,j))*dt*fdp1)
	    
   end if

    if (fdm1>=1d0-fc.and.fd>=fc.and.fdp1<=fc) then
        flux_y(i,j)=min(fd*dy,max(0d0,abs(v(i,j))*dt-(1-fd)*dy))
		
   end if

    if (fdm1>=fc.and.fd>=fc.and.fdp1>=fc) then
        flux_y(i,j)=abs(v(i,j))*dt*fd
		
   end if

    if (fdm1>=1-fc.and.fd>=fc.and.fdp1>=fc) then
        flux_y(i,j)=min(fd*dy,fdp1*abs(v(i,j))*dt+max((1-fdp1)*abs(v(i,j))*dt-(1d0-fd)*dy,0d0))
		
   end if

    if (fdm1>=fc.and.fd>=fc.and.fdp1>=1-fc) then
       flux_y(i,j)=min(fdp1*abs(v(i,j))*dt, fd*dy-fdm1*(dy-abs(v(i,j))*dt))
        	
   end if

   us=(fdm1<fc.and.fd<fc.and.fdp1<fc).or.(fdm1>fc.and.fd>fc.and.fdp1>fc).or.&
	   (fdm1<fc.and.fd<fc).or.(fd<fc.and.fdp1<fc)
	if(us) then
	   flux_y(i,j)=abs(v(i,j))*dt*fd
	end if
end do 
end do
do i=2,L1-1
 do j=2,M1-1   
   fs(i,j)=fsx(i,j)*(1d0+dt/dy*(v(i,j)-v(i,j-1)))+1d0/dy*(sign2(v(i,j-1))*flux_y(i,j-1)-sign2(v(i,j))*flux_y(i,j)) !!!!!!!!
   
end do
end do
else
   

   !solving fluxes for y-axis
do i=2,L1-1
 do j=2,M1-1
   if(v(i,j)>0d0) then !!!!!!!!!!!!!!!!
     fd=fs(i,j)
	 fdm1=fs(i,j-1)
	 fdp1=fs(i,j+1)
	 
   end if

   if(v(i,j)<0d0) then
     fd=fs(i,j)
	 fdm1=fs(i,j+1)
	 fdp1=fs(i,j-1)
	
   end if

   !if(v(i,j)==0d0) exit

   if (fdm1<=fc.and.fd>=fc.and.fdp1>=fc) then
       flux_y(i,j)=min(fd*dy,abs(v(i,j))*dt*fdp1)
       
   end if

    if (fdm1>=1d0-fc.and.fd>=fc.and.fdp1<=fc) then
        flux_y(i,j)=min(fd*dy,max(0d0,abs(v(i,j))*dt-(1-fd)*dy))
        
   end if

    if (fdm1>=fc.and.fd>=fc.and.fdp1>=fc) then
        flux_y(i,j)=abs(v(i,j))*dt*fd
        
   end if

    if (fdm1>=1-fc.and.fd>=fc.and.fdp1>=fc) then
        flux_y(i,j)=min(fd*dy,abs(v(i,j))*dt*fdp1+max((1-fdp1)*abs(v(i,j))*dt-(1d0-fd)*dy,0d0))
		
   end if

    if (fdm1>=fc.and.fd>=fc.and.fdp1>=1-fc) then
       flux_y(i,j)=min(fdp1*abs(v(i,j))*dt, fd*dy-fdm1*(dy-abs(v(i,j))*dt))
       
   end if

   us=(fdm1<fc.and.fd<fc.and.fdp1<fc).or.(fdm1>fc.and.fd>fc.and.fdp1>fc).or.&
	   (fdm1<fc.and.fd<fc).or.(fd<fc.and.fdp1<fc)
	if(us) then
	   flux_y(i,j)=abs(v(i,j))*dt*fd
	end if
end do
end do

do i=2,L1-1
 do j=2,M1-1
   fsy(i,j)=(fs(i,j)+1d0/dy*(sign2(v(i,j-1))*flux_y(i,j-1)-sign2(v(i,j))*flux_y(i,j)))/(1d0-dt/dy*(v(i,j)-v(i,j-1)))
   
end do
end do
      !solving fluxes for x-axis

do i=2,L1-1
 do j=2,M1-1
   if(u(i,j)>0d0) then
     fd=fsy(i,j)
	 fdm1=fsy(i-1,j)
	 fdp1=fsy(i+1,j)
     
   end if

   if(u(i,j)<0d0) then
     fd=fsy(i,j)
	 fdm1=fsy(i+1,j)
	 fdp1=fsy(i-1,j)
	
   end if
   
   !if(u(i,j)==0d0) exit

   if (fdm1<=fc.and.fd>=fc.and.fdp1>=fc) then
       flux_x(i,j)=min(fd*dx,abs(u(i,j))*dt*fdp1)
       
       
   end if

    if (fdm1>=1d0-fc.and.fd>=fc.and.fdp1<=fc) then
        flux_x(i,j)=min(fd*dx,max(0d0,abs(u(i,j))*dt-(1-fd)*dx))
		
    end if

    if (fdm1>=fc.and.fd>=fc.and.fdp1>=fc) then
        flux_x(i,j)=abs(u(i,j))*dt*fd  !!!!!!!error
		
    end if

    if (fdm1>=1-fc.and.fd>=fc.and.fdp1>=fc) then
        flux_x(i,j)=min(fd*dx,abs(u(i,j))*dt*fdp1+max((1-fdp1)*abs(u(i,j))*dt-(1d0-fd)*dx,0d0))
		
    end if

    if (fdm1>=fc.and.fd>=fc.and.fdp1>=1-fc) then
       flux_x(i,j)=min(fdp1*abs(u(i,j))*dt, fd*dx-fdm1*(dx-abs(u(i,j))*dt))
	   
    end if

	us=(fdm1<fc.and.fd<fc.and.fdp1<fc).or.(fdm1>fc.and.fd>fc.and.fdp1>fc).or.&
	   (fdm1<fc.and.fd<fc).or.(fd<fc.and.fdp1<fc)
	if(us) then
	   flux_x(i,j)=abs(u(i,j))*dt*fd
	end if
 end do
end do

do i=2,L1-1
 do j=2,M1-1
   fs(i,j)=fsy(i,j)*(1d0+dt/dx*(u(i,j)-u(i-1,j)))+1d0/dx*(sign2(u(i-1,j))*flux_x(i-1,j)-sign2(u(i,j))*flux_x(i,j))
   

 end do
end do
end if
end