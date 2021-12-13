subroutine advect_da2(u,v,fs,fb,mark,ax,ay,snormx,snormy,c,x,y,L1,M1,dx,dy,dt,t,itimestep)
!plic donor acceptor method y-x sweep not working
use vof2d
implicit none
integer,intent(in)::L1,M1,itimestep
integer,intent(inout)::mark(ni,nj)
real(8),intent(in)::u(ni,nj),v(ni,nj),ax(ni,nj),ay(ni,nj),fb(ni,nj),x(ni),y(nj)
real(8),intent(in)::t,dt,dx,dy
real(8),intent(inout)::fs(ni,nj)
real(8),intent(inout):: snormx(ni,nj),snormy(ni,nj),c(ni,nj)
real(8),external::sign2
real(8),parameter::fc=0.01d0,fss=0.4d0,sc=0.5d0
integer i,j,fptr
real(8) fa,fd,fad,cnw,csw,cne,cse,sp,xe_new,xw_new,fba,fbd,fbad,cf,fdd,fbdd
!real(8) snormx(ni,nj),snormy(ni,nj),c(ni,nj)
real(8) flux_x(ni,nj),flux_y(ni,nj),fsx(ni,nj),fsy(ni,nj),slope(ni,nj)
real(8) xx(ni,nj),yy(ni,nj),dxy(ni,nj),dyx(ni,nj),dyxa,dxya,dyxd,dxyd
character*200 filename1,szNumber
logical alpha
fptr=125
Write(szNumber,'(i6.6)') itimestep
filename1=Trim("../voffunc2/test/")//trim('plic')//trim(szNumber)//trim('.txt')
!open( fptr,file=filename1,status='new')	  
alpha=.true.

flux_x=0d0; flux_y=0d0; fsx=0d0; fsy=0d0;c=0d0;slope=0d0
if(alpha)then
!if(mod(itimestep,2)==0)then
write(fptr,'(A)') "solving tr.equation x-y sweep,y-x sweep"
write(fptr,'(A)') "x-y sweep"
call PRNT(fs,"fs_old",fptr,L1,M1)
write(fptr,*) "   "
!x-y sweep
  call normal(.true.,L1,M1,x,y,fs,mark,snormx,snormy)
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fs(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		!if(abs(snormx(i,j))==1d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=huge(1d0)
	    !else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=0d0 !slope(i,j)=huge(1d0)        
	    !else if(fb(i,j)==fs(i,j)) then
		!  slope(i,j)=huge(1d0)
		!else 
	    !  slope(i,j)=-(snormx(i,j)/snormy(i,j))
	    !end if
	  end if
	end do
  end do
  !solving fluxes for x-y axis
do i=2,L1-1
  do j=2,M1-1
   xx(i,j)=dx*(fs(i-1,j)+fs(i,j)+fs(i+1,j))
   yy(i,j)=dy*(fs(i,j-1)+fs(i,j)+fs(i,j+1))
  end do
end do

do i=2,L1-1
  do j=2,M1-1
   dxy(i,j)=2d0*(xx(i,j+1)-xx(i,j-1))/(4d0*dy)
   dyx(i,j)=2d0*(xx(i+1,j)-xx(i-1,j))/(4d0*dy)
  end do
end do


do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i-1,j)/=B.or.mark(i+1,j)/=B) then 
	    if(u(i,j)>=0d0) then
          fd=fs(i,j)	 
	      fa=fs(i+1,j)
		  fdd=fs(i-1,j)
          fbd=fb(i,j)
		  fba=fb(i+1,j)
		  fbdd=fb(i-1,j)
		  dyxa=dyx(i+1,j)
		  dxya=dxy(i+1,j)
		  dyxd=dyx(i,j)
		  dxyd=dxy(i,j)
		  
        end if
        if(u(i,j)<=0d0) then
          fd=fs(i+1,j)	 
	      fa=fs(i,j)
		  fdd=fs(i+2,j)
		  fbd=fb(i+1,j)	 
	      fba=fb(i,j)
          fbdd=fb(i+2,j)
          dyxa=dyx(i,j)
		  dxya=dxy(i,j)
		  dyxd=dyx(i+1,j)
		  dxyd=dxy(i+1,j)

        end if

		if(abs(dyxd)<=abs(dxyd).and.fa/=0d0.and.fdd/=0d0) then
		   fad=fd
		   fbad=fbd
		else if(abs(dyxa)>abs(dxya).or.fa==0d0.or.fdd==0d0) then
           fad=fa
		   fbad=fba 
		end if
        
        
		
        if(u(i,j)==0d0) then
           flux_x(i,j)=0d0
		else
		  cf=max((fbad-fad)/fbad*abs(ax(i,j)*u(i,j))*dt-(fbd-fd)*dx,0d0)
          flux_x(i,j)=min(fad/fbad*abs(ax(i,j)*u(i,j))*dt+cf,fd*dx)
		end if 
	 end if
   end do
 end do



  
 !solving transport equation for x-sweep
 do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    fsx(i,j)=(fs(i,j)+1d0/dx*(sign2(u(i-1,j))*flux_x(i-1,j)&
		          -sign2(u(i,j))*flux_x(i,j)))/(1d0-dt/dx*(u(i,j)-u(i-1,j)))
	  end if
	end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(fptr,'(A)') "before x-y sweep"
call PRNT(snormx,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(snormy,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(c,"const",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(slope,"slope",fptr,L1,M1)
write(fptr,*) "   "
write(fptr,'(A)') "solving x-y sweep"
call PRNT(flux_x,"flux_x",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(fsx,"fsx",fptr,L1,M1)
write(fptr,*) "   "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !y-x sweep
  write(fptr,'(A)') "y-x sweep"
  slope=0d0
  call marking(fb,fsx,mark,L1,M1)
  call normal(.true.,L1,M1,x,y,fsx,mark,snormx,snormy)
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fsx(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fsx(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fsx(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		!if(abs(snormx(i,j))==1d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=huge(1d0)
	    !else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=0d0 !slope(i,j)=huge(1d0)! 
	    !else if(fb(i,j)==fs(i,j)) then
		!  slope(i,j)=huge(1d0)
		!else
	    !  slope(i,j)=-(snormx(i,j)/snormy(i,j))
	    !end if
        
	  end if
	end do
  end do

 do i=2,L1-1
  do j=2,M1-1
   xx(i,j)=dx*(fsx(i-1,j)+fsx(i,j)+fsx(i+1,j))
   yy(i,j)=dy*(fsx(i,j-1)+fsx(i,j)+fsx(i,j+1))
  end do
end do

do i=2,L1-1
  do j=2,M1-1
   dxy(i,j)=2d0*(xx(i,j+1)-xx(i,j-1))/(4d0*dy)
   dyx(i,j)=2d0*(xx(i+1,j)-xx(i-1,j))/(4d0*dy)
  end do
end do



  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i,j-1)/=B.or.mark(i,j+1)/=B) then
	    if(v(i,j)>0d0) then
          fd=fsx(i,j)	 
	      fa=fsx(i,j+1)
		  fdd=fsx(i,j-1)
		  fbd=fb(i,j)	 
	      fba=fb(i,j+1)
		  fbdd=fb(i,j-1)
		  dyxa=dyx(i,j+1)
		  dxya=dxy(i,j+1)
		  dyxd=dyx(i,j)
		  dxyd=dxy(i,j)	     
        end if
        if(v(i,j)<0d0) then
          fd=fsx(i,j+1)	 
	      fa=fsx(i,j)
		  fdd=fsx(i,j+2)
		  fbd=fb(i,j+1)	 
	      fba=fb(i,j)
		  fbdd=fb(i,j+2)
		  dyxa=dyx(i,j)
		  dxya=dxy(i,j)
		  dyxd=dyx(i,j+1)
		  dxyd=dxy(i,j+1)	      
        end if

		        

		if(abs(dyxd)>=abs(dxyd).and.fa/=0d0.and.fdd/=0d0) then
		   fad=fd
		   fbad=fbd
		else if(abs(dyxa)<abs(dxya).or.fa==0d0.or.fdd==0d0) then
           fad=fa
		   fbad=fba 
		end if

        if(v(i,j)==0d0) then
		  flux_y(i,j)=0d0
		else
		  cf=max((fbad-fad)/fbad*abs(ay(i,j)*v(i,j))*dt-(fbd-fd)*dy,0d0)
          flux_y(i,j)=min(fad/fbad*abs(ay(i,j)*v(i,j))*dt+cf,fd*dy)
	    end if
	 end if
   end do
 end do
   
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    fs(i,j)=fsx(i,j)*(1d0+dt/dy*(v(i,j)-v(i,j-1))) &
		        +1d0/dy*(sign2(v(i,j-1))*flux_y(i,j-1)-sign2(v(i,j))*flux_y(i,j))
	  end if
   end do
 end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(fptr,'(A)') "before y-x sweep"
call PRNT(snormx,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(snormy,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(c,"const",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(slope,"slope",fptr,L1,M1)
write(fptr,*) "   "
write(fptr,'(A)') "solving y-x sweep"
call PRNT(flux_y,"flux_y",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(fs,"fs",fptr,L1,M1)
write(fptr,*) "   "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  slope=0d0
  call marking(fb,fs,mark,L1,M1)
  call normal(.true.,L1,M1,x,y,fs,mark,snormx,snormy)
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fs(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		!if(abs(snormx(i,j))==1d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=huge(1d0)
	    !else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=0d0
	    !else
	    !  slope(i,j)=-(snormx(i,j)/snormy(i,j))
	    !end if
	  end if
	end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(fptr,'(A)') "final PLIC reconstruction after y-x sweep"
call PRNT(snormx,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(snormy,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(c,"const",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(slope,"slope",fptr,L1,M1)
write(fptr,*) "   "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
else

!y-x sweep
write(fptr,'(A)') "solving tr.equation y-x sweep,x-y sweep"
write(fptr,'(A)') "y-x sweep"
call PRNT(fs,"fs_old",fptr,L1,M1)
write(fptr,*) "   "
!call marking(fb,fs,mark,L1,M1)
  call normal(.true.,L1,M1,x,y,fs,mark,snormx,snormy)
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fs(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		!if(abs(snormx(i,j))==1d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=huge(1d0)
	    !else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=0d0
	    !else
	     ! slope(i,j)=-(snormx(i,j)/snormy(i,j))
	    !end if
	  end if
	end do
  end do


 do i=2,L1-1
  do j=2,M1-1
   xx(i,j)=dx*(fs(i-1,j)+fs(i,j)+fs(i+1,j))
   yy(i,j)=dy*(fs(i,j-1)+fs(i,j)+fs(i,j+1))
  end do
end do

do i=2,L1-1
  do j=2,M1-1
   dxy(i,j)=2d0*(xx(i,j+1)-xx(i,j-1))/(4d0*dy)
   dyx(i,j)=2d0*(xx(i+1,j)-xx(i-1,j))/(4d0*dy)
  end do
end do



  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i,j-1)/=B.or.mark(i,j+1)/=B) then
	    if(v(i,j)>0d0) then
          fd=fsx(i,j)	 
	      fa=fsx(i,j+1)
		  fdd=fsx(i,j-1)
		  fbd=fb(i,j)	 
	      fba=fb(i,j+1)
		  fbdd=fb(i,j-1)
		  dyxa=dyx(i,j+1)
		  dxya=dxy(i,j+1)
		  dyxd=dyx(i,j)
		  dxyd=dxy(i,j)		  
        end if
        if(v(i,j)<0d0) then
          fd=fsx(i,j+1)	 
	      fa=fsx(i,j)
		  fdd=fsx(i,j+2)
		  fbd=fb(i,j+1)	 
	      fba=fb(i,j)
		  fbdd=fb(i,j+2)
		  dyxa=dyx(i,j)
		  dxya=dxy(i,j)
		  dyxd=dyx(i,j+1)
		  dxyd=dxy(i,j+1)	      
        end if

       
        

        if(abs(dyxd)<=abs(dxyd).and.fa/=0d0.and.fdd/=0d0) then
		   fad=fd
		   fbad=fbd
		else if(abs(dyxa)>abs(dxya).or.fa==0d0.or.fdd==0d0) then
           fad=fa
		   fbad=fba 
		end if

        if(v(i,j)==0d0) then
		  flux_y(i,j)=0d0
		else
          cf=max((fbad-fad)/fbad*abs(ay(i,j)*v(i,j))*dt-(fbd-fd)*dy,0d0)
          flux_y(i,j)=min(fad/fbad*abs(ay(i,j)*v(i,j))*dt+cf,fd*dy)
	    end if

	 end if
   end do
 end do
   
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    fsy(i,j)=(fs(i,j)+1d0/dy*(sign2(v(i,j-1))*flux_y(i,j-1) &
		          -sign2(v(i,j))*flux_y(i,j)))/(1d0-dt/dy*(v(i,j)-v(i,j-1)))
	  end if
   end do
 end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(fptr,'(A)') "before y-x sweep"
call PRNT(snormx,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(snormy,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(c,"const",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(slope,"slope",fptr,L1,M1)
write(fptr,*) "   "
write(fptr,'(A)') "solving y-x sweep"
call PRNT(flux_y,"flux_y",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(fsy,"fsy",fptr,L1,M1)
write(fptr,*) "   "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!solving flux for x-axis
slope=0d0
 call marking(fb,fsy,mark,L1,M1)
 call normal(.true.,L1,M1,x,y,fsy,mark,snormx,snormy)
 do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fsy(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fsy(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fsy(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		!if(abs(snormx(i,j))==1d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=huge(1d0)
	    !else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=0d0
	    !else
	    !  slope(i,j)=-(snormx(i,j)/snormy(i,j))
	    !end if
	  end if
	end do
  end do

 do i=2,L1-1
  do j=2,M1-1
   xx(i,j)=dx*(fsy(i-1,j)+fsy(i,j)+fsy(i+1,j))
   yy(i,j)=dy*(fsy(i,j-1)+fsy(i,j)+fsy(i,j+1))
  end do
end do

do i=2,L1-1
  do j=2,M1-1
   dxy(i,j)=2d0*(xx(i,j+1)-xx(i,j-1))/(4d0*dy)
   dyx(i,j)=2d0*(xx(i+1,j)-xx(i-1,j))/(4d0*dy)
  end do
end do


  !solving fluxes for x-y axis
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i-1,j)/=B.or.mark(i+1,j)/=B) then
	    if(u(i,j)>=0d0) then
          fd=fsy(i,j)	 
	      fa=fsy(i+1,j)
		  fdd=fsy(i-1,j)
          fbd=fb(i,j)
		  fba=fb(i+1,j)
		  fbdd=fb(i-1,j)
		  dyxa=dyx(i+1,j)
		  dxya=dxy(i+1,j)
		  dyxd=dyx(i,j)
		  dxyd=dxy(i,j)		  
        end if
        if(u(i,j)<=0d0) then
          fd=fsy(i+1,j)	 
	      fa=fsy(i,j)
		  fdd=fsy(i+2,j)
		  fbd=fb(i+1,j)	 
	      fba=fb(i,j)
          fbdd=fb(i+2,j)
		  dyxa=dyx(i,j)
		  dxya=dxy(i,j)
		  dyxd=dyx(i+1,j)
		  dxyd=dxy(i+1,j)
        end if

        if(abs(dyxd)>=abs(dxyd).and.fa/=0d0.and.fdd/=0d0) then
		   fad=fd
		   fbad=fbd
		else if(abs(dyxa)<abs(dxya).or.fa==0d0.or.fdd==0d0) then
           fad=fa
		   fbad=fba 
		end if

       
        if(u(i,j)==0d0) then
           flux_x(i,j)=0d0
		else
          cf=max((fbad-fad)/fbad*abs(ax(i,j)*u(i,j))*dt-(fbd-fd)*dx,0d0)
          flux_x(i,j)=min(fad/fbad*abs(ax(i,j)*u(i,j))*dt+cf,fd*dx)
		end if 
	 end if
   end do
 end do
 !solving transport equation for x-sweep
 do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    fs(i,j)=fsy(i,j)*(1d0+dt/dx*(u(i,j)-u(i-1,j)))&
		        +1d0/dx*(sign2(u(i-1,j))*flux_x(i-1,j)-sign2(u(i,j))*flux_x(i,j))
	  end if
	end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(fptr,'(A)') "before x-y sweep"
call PRNT(snormx,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(snormy,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(c,"const",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(slope,"slope",fptr,L1,M1)
write(fptr,*) "   "
write(fptr,'(A)') "solving x-y sweep"
call PRNT(flux_x,"flux_x",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(fs,"fs",fptr,L1,M1)
write(fptr,*) "   "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  slope=0d0
  call marking(fb,fs,mark,L1,M1)
  call normal(.true.,L1,M1,x,y,fs,mark,snormx,snormy)
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fs(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		!if(abs(snormx(i,j))==1d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=huge(1d0)
	    !else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
        !  slope(i,j)=0d0
	    !else
	    !  slope(i,j)=-(snormx(i,j)/snormy(i,j))
	    !end if
	  end if
	end do
  end do
write(fptr,'(A)') "final PLIC reconstruction after x-y sweep"
call PRNT(snormx,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(snormy,"nx",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(c,"const",fptr,L1,M1)
write(fptr,*) "   "
call PRNT(slope,"slope",fptr,L1,M1)
write(fptr,*) "   "
end if

close(fptr)
end subroutine advect_da2

