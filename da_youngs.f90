subroutine advect_da3(u,v,fs,fb,mark,ax,ay,snormx,snormy,c,x,y,L1,M1,dx,dy,dt,t,itimestep)
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
real(8) fa,fd,fad,cnw,csw,cne,cse,sp,xe_new,xw_new,fba,fbd,fbad,cf
!real(8) snormx(ni,nj),snormy(ni,nj),c(ni,nj)
real(8) flux_x(ni,nj),flux_y(ni,nj),fsx(ni,nj),fsy(ni,nj),slope(ni,nj)
real(8) cp,nx,ny,xw,xe,ys,yn,vol
character*200 filename1,szNumber
logical alpha
!fptr=125
!Write(szNumber,'(i6.6)') itimestep
!filename1=Trim("../voffunc2/test/")//trim('plic')//trim(szNumber)//trim('.txt')
!open( fptr,file=filename1,status='new')	  
alpha=.true.

flux_x=0d0; flux_y=0d0; fsx=0d0; fsy=0d0;c=0d0;slope=0d0
if(alpha)then
!if(mod(itimestep,2)==0)then
!write(fptr,'(A)') "solving tr.equation x-y sweep,y-x sweep"
!write(fptr,'(A)') "x-y sweep"
!call PRNT(fs,"fs_old",fptr,L1,M1)
!write(fptr,*) "   "
!x-y sweep
  call normal(.true.,L1,M1,x,y,fs,mark,snormx,snormy)
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    if(mark(i,j)==S) then 
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fs(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		end if

		if(abs(snormx(i,j))==1d0.and.snormy(i,j)==0d0) then
          slope(i,j)=huge(1d0)
	    else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
          slope(i,j)=0d0 
		else 
	      slope(i,j)=-(snormx(i,j)/snormy(i,j))
	    end if
	  end if
	end do
  end do
  !solving fluxes for x-y axis

  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i-1,j)/=B.or.mark(i+1,j)/=B) then 
	    if(u(i,j)>0d0) then
          fd=fs(i,j)	 
	      fa=fs(i+1,j)
          fbd=fb(i,j)
		  fba=fb(i+1,j)
		  
        end if
        if(u(i,j)<0d0) then
          fd=fs(i+1,j)	 
	      fa=fs(i,j)
		  fbd=fb(i+1,j)	 
	      fba=fb(i,j)

        end if

		if(u(i,j)>0d0) then
          sp=slope(i,j)
		  nx=snormx(i,j)
		  ny=snormy(i,j)
		  cp=c(i,j)
		  xe=x(i)
		  xw=x(i)-u(i,j)*dt
		  yn=y(j)
		  ys=y(j-1)
		end if

        if(u(i,j)<0d0) then
          sp=slope(i+1,j)
          nx=snormx(i+1,j)
		  ny=snormy(i+1,j)
		  cp=c(i+1,j)
		  xw=x(i)
		  xe=x(i)+abs(u(i,j))*dt
		  yn=y(j)
		  ys=y(j-1)
		end if

        if(abs(sp)>=sc) then   
          fad=fa
		  fbad=fba
        else
          fad=fd
		  fbad=fbd
        end if
		
         if(u(i,j)==0d0) then
           flux_x(i,j)=0d0
		else
		  if(fad<=0d0+1d-5) then
		    vol=0d0
		  else if(fad>=1d0-1d-5)then
		     vol=abs(u(i,j))*dt*dy
		  else
		    call plc(nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse)
		    call value2(cp,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,vol)
			
		  end if
		  cf=max((fbad-fad)/fbad*ax(i,j)*vol/dy-(fbd-fd)*dx,0d0)
          flux_x(i,j)=min(fad/fbad*ax(i,j)*vol/dy+cf,fd*dx)
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
!write(fptr,'(A)') "before x-y sweep"
!call PRNT(snormx,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(snormy,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(c,"const",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(slope,"slope",fptr,L1,M1)
!write(fptr,*) "   "
!write(fptr,'(A)') "solving x-y sweep"
!call PRNT(flux_x,"flux_x",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(fsx,"fsx",fptr,L1,M1)
!write(fptr,*) "   "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !y-x sweep
!  write(fptr,'(A)') "y-x sweep"
  slope=0d0
  call marking(fb,fsx,mark,L1,M1)
  call normal(.true.,L1,M1,x,y,fsx,mark,snormx,snormy)
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    if(mark(i,j)==S) then 
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fsx(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fsx(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fsx(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		end if
		if(abs(snormy(i,j))==1d0.and.snormx(i,j)==0d0) then
          slope(i,j)=huge(1d0)
	    else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
          slope(i,j)=0d0 
		else
	      slope(i,j)=-(snormy(i,j)/snormx(i,j))
	    end if
        
	  end if
	end do
  end do

  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i,j-1)/=B.or.mark(i,j+1)/=B) then
	    if(v(i,j)>0d0) then
          fd=fsx(i,j)	 
	      fa=fsx(i,j+1)
		  fbd=fb(i,j)	 
	      fba=fb(i,j+1)	     
        end if
        if(v(i,j)<0d0) then
          fd=fsx(i,j+1)	 
	      fa=fsx(i,j)
		  fbd=fb(i,j+1)	 
	      fba=fb(i,j)	      
        end if

		if(v(i,j)>0d0) then
          sp=slope(i,j)
		  nx=snormx(i,j)
		  ny=snormy(i,j)
		  cp=c(i,j)
		  xe=x(i)
		  xw=x(i-1)
		  yn=y(j)
		  ys=y(j)-v(i,j)*dt
		  
		end if

		if(v(i,j)<0d0) then          
		  sp=slope(i,j+1)
		  nx=snormx(i,j+1)
		  ny=snormy(i,j+1)
		  cp=c(i,j+1)
		  xe=x(i)
		  xw=x(i-1)
		  ys=y(j)
		  yn=y(j)+abs(v(i,j))*dt
		end if        

		if(abs(sp)>=sc) then   
          fad=fa
		  fbad=fba
        else
          fad=fd
		  fbad=fbd
        end if

        if(v(i,j)==0d0) then
		  flux_y(i,j)=0d0
		else
		  if(fad<=0d0+1d-5) then
		    vol=0d0
		  else if(fad>=1d0-1d-5)then
		    vol=abs(v(i,j))*dt*dx
		  else 
		   call plc(nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse)
		   call value2(cp,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,vol)
		   cf=max((fbad-fad)/fbad*ay(i,j)*vol/dx-(fbd-fd)*dy,0d0)
           flux_y(i,j)=min(fad/fbad*ay(i,j)*vol/dx+cf,fd*dy)
		 end if
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
!write(fptr,'(A)') "before y-x sweep"
!call PRNT(snormx,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(snormy,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(c,"const",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(slope,"slope",fptr,L1,M1)
!write(fptr,*) "   "
!write(fptr,'(A)') "solving y-x sweep"
!call PRNT(flux_y,"flux_y",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(fs,"fs",fptr,L1,M1)
!write(fptr,*) "   "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  slope=0d0
  call marking(fb,fs,mark,L1,M1)
  call normal(.true.,L1,M1,x,y,fs,mark,snormx,snormy)
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    if(mark(i,j)==S) then 
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fs(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		end if
		if(abs(snormx(i,j))==1d0.and.snormy(i,j)==0d0) then
          slope(i,j)=huge(1d0)
	    else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
          slope(i,j)=0d0
	    else
	      slope(i,j)=-(snormx(i,j)/snormy(i,j))
	    end if
	  end if
	end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(fptr,'(A)') "final PLIC reconstruction after y-x sweep"
!call PRNT(snormx,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(snormy,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(c,"const",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(slope,"slope",fptr,L1,M1)
!write(fptr,*) "   "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
else

!y-x sweep
!write(fptr,'(A)') "solving tr.equation y-x sweep,x-y sweep"
!write(fptr,'(A)') "y-x sweep"
!call PRNT(fs,"fs_old",fptr,L1,M1)
!write(fptr,*) "   "
  call marking(fb,fs,mark,L1,M1)
  call normal(.true.,L1,M1,x,y,fs,mark,snormx,snormy)
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    if(mark(i,j)==S) then 
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fs(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		end if
		if(abs(snormy(i,j))==1d0.and.snormx(i,j)==0d0) then
          slope(i,j)=huge(1d0)
	    else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
          slope(i,j)=0d0 
		else
	      slope(i,j)=-(snormy(i,j)/snormx(i,j))
	    end if
	  end if
	end do
  end do

  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i,j-1)/=B.or.mark(i,j+1)/=B) then
	    if(v(i,j)>0d0) then
          fd=fs(i,j)	 
	      fa=fs(i,j+1)
		  fbd=fb(i,j)	 
	      fba=fb(i,j+1)	     
        end if
        if(v(i,j)<0d0) then
          fbd=fb(i,j+1)	 
	      fba=fb(i,j)	      
        end if

        if(v(i,j)>0d0) then
          sp=slope(i,j)
          nx=snormx(i,j)
		  ny=snormy(i,j)
		  cp=c(i,j)
		  xe=x(i)
		  xw=x(i-1)
		  yn=y(j)
		  ys=y(j)-v(i,j)*dt
		end if

		if(v(i,j)<0d0) then
          sp=slope(i,j+1)
		  nx=snormx(i,j+1)
		  ny=snormy(i,j+1)
		  cp=c(i,j+1)
		  xe=x(i)
		  xw=x(i-1)
		  ys=y(j)
		  yn=y(j)+abs(v(i,j))*dt
		end if

        

        if(abs(sp)>=sc) then   
          fad=fa
		  fbad=fba
        else
          fad=fd
		  fbad=fbd
        end if

        if(v(i,j)==0d0) then
		  flux_y(i,j)=0d0
		else
		  if(fad<=0d0+1d-5) then
		    vol=0d0
		  else if(fad>=1d0-1d-5)then
		     vol=abs(v(i,j))*dt*dx
		  else 
           call plc(nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse)
		   call value2(cp,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,vol)
		   cf=max((fbad-fad)/fbad*ay(i,j)*vol/dx-(fbd-fd)*dy,0d0)
           flux_y(i,j)=min(fad/fbad*ay(i,j)*vol/dx+cf,fd*dy)
		  end if
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
!write(fptr,'(A)') "before y-x sweep"
!call PRNT(snormx,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(snormy,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(c,"const",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(slope,"slope",fptr,L1,M1)
!write(fptr,*) "   "
!write(fptr,'(A)') "solving y-x sweep"
!call PRNT(flux_y,"flux_y",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(fsy,"fsy",fptr,L1,M1)
!write(fptr,*) "   "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!solving flux for x-axis
slope=0d0;snormx=0d0;snormy=0d0;c=0d0
 call marking(fb,fsy,mark,L1,M1)
 call normal(.true.,L1,M1,x,y,fsy,mark,snormx,snormy)
   do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    if(mark(i,j)==S) then 
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fsy(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fsy(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fsy(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		end if
		if(abs(snormx(i,j))==1d0.and.snormy(i,j)==0d0) then
          slope(i,j)=huge(1d0)
	    else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
          slope(i,j)=0d0
	    else
	      slope(i,j)=-(snormx(i,j)/snormy(i,j))
	    end if
	  end if
	end do
  end do !!!!gluck v etom cycle
  !solving fluxes for x-y axis
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i-1,j)/=B.or.mark(i+1,j)/=B) then
	    if(u(i,j)>0d0) then
          fd=fsy(i,j)	 
	      fa=fsy(i+1,j)
		  fbd=fb(i,j)	 
	      fba=fb(i+1,j)	     
        end if
        if(u(i,j)<0d0) then
          fd=fsy(i+1,j)	 
	      fa=fsy(i,j)
		  fbd=fb(i+1,j)	 
	      fba=fb(i,j)	      
        end if

        if(u(i,j)>0d0) then
          sp=slope(i,j)
		  nx=snormx(i,j)
		  ny=snormy(i,j)
		  cp=c(i,j)
		  xe=x(i)
		  xw=x(i)-u(i,j)*dt
		  yn=y(j)
		  ys=y(j-1)
		end if

        if(u(i,j)<0d0) then
          sp=slope(i+1,j)
          sp=slope(i+1,j)
          nx=snormx(i+1,j)
		  ny=snormy(i+1,j)
		  cp=c(i+1,j)
		  xw=x(i)
		  xe=x(i)+abs(u(i,j))*dt
		  yn=y(j)
		  ys=y(j-1)
		end if

        if(abs(sp)>=sc) then   
          fad=fa
		  fbad=fba
        else
		  fad=fd
          fbad=fbd
        end if
        if(u(i,j)==0d0) then
           flux_x(i,j)=0d0
		else
		  if(fad<=0d0+1d-5) then
		    vol=0d0
		   else if(fad>=1d0-1d-5)then
		    vol=abs(u(i,j))*dt*dy
		  else 
		    call plc(nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse)
		    call value2(cp,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,vol)
            cf=max((fbad-fad)/fbad*ax(i,j)*vol/dy-(fbd-fd)*dx,0d0)
            flux_x(i,j)=min(fad/fbad*ax(i,j)*vol/dy+cf,fd*dx)
		  end if
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
!write(fptr,'(A)') "before x-y sweep"
!call PRNT(snormx,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(snormy,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(c,"const",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(slope,"slope",fptr,L1,M1)
!write(fptr,*) "   "
!write(fptr,'(A)') "solving x-y sweep"
!call PRNT(flux_x,"flux_x",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(fs,"fs",fptr,L1,M1)
!write(fptr,*) "   "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  slope=0d0
  call marking(fb,fs,mark,L1,M1)
  call normal(.true.,L1,M1,x,y,fs,mark,snormx,snormy)
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    if(mark(i,j)==S) then 
		if(fb(i,j)<1d0) then
		  call reset_xy(x(i-1),x(j),y(j),y(j-1),fs(i,j),fb(i,j),ax(i-1,j),&
		                 ax(i,j),ay(i,j),ay(i,j-1),xe_new,xw_new)
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),xw_new,xe_new,y(j),y(j-1),c(i,j))
		else
		  call planec(fs(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),c(i,j))
		end if
		end if
		if(abs(snormx(i,j))==1d0.and.snormy(i,j)==0d0) then
          slope(i,j)=huge(1d0)
	    else if(snormx(i,j)==0d0.and.snormy(i,j)==0d0) then
          slope(i,j)=0d0
	    else
	      slope(i,j)=-(snormx(i,j)/snormy(i,j))
	    end if
	  end if
	end do
  end do
!write(fptr,'(A)') "final PLIC reconstruction after x-y sweep"
!call PRNT(snormx,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(snormy,"nx",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(c,"const",fptr,L1,M1)
!write(fptr,*) "   "
!call PRNT(slope,"slope",fptr,L1,M1)
!write(fptr,*) "   "
end if

!close(fptr)
end subroutine advect_da3
