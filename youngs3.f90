subroutine advect_da5(u,v,fs,fb,mark,ax,ay,snormx,snormy,c,x,y,L1,M1,dx,dy,dt,t,itimestep)
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
real(8) flux_x(ni,nj),flux_y(ni,nj),fsx(ni,nj),fsy(ni,nj),slope(ni,nj),vn,vofsum
real(8) cp,nx,ny,xw,xe,ys,yn,va,vd,vad,nxa,nxd,nya,nyd,cpd,cpa,vfd,vfa,vfad,vfl(ni,nj),vfr(ni,nj),vol(ni,nj)      !v1,v2,vf1,vf2
character*200 filename1,szNumber
logical alpha
!fptr=125
!Write(szNumber,'(i6.6)') itimestep
!filename1=Trim("../voffunc2/test/")//trim('plic')//trim(szNumber)//trim('.txt')
!open( fptr,file=filename1,status='new')	  
alpha=.true.

flux_x=0d0; flux_y=0d0; fsx=0d0; fsy=0d0;c=0d0;slope=0d0
!if(alpha)then
if(mod(itimestep,2)==0)then
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

      if(fs(i,j)<=0d0+1d-5) then	    
		vfl(i,j)=0d0
		vfr(i,j)=0d0
		vol(i,j)=0d0
	  else if(fs(i,j)>=1d0-1d-5)then
	    vfr(i,j)=abs(u(i,j))*dt*dy
		vfl(i,j)=abs(u(i,j))*dt*dy
		vol(i,j)=dx*dy
	  else
	    if(mark(i,j)/=E.or.mark(i,j)/=F)then
		    call plc(snormx(i,j),snormy(i,j),x(i)-abs(u(i,j))*dt,x(i),y(j),y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i)-abs(u(i,j))*dt,x(i),y(j),y(j-1),cnw,csw,cne,cse,vfr(i,j))

			call plc(snormx(i,j),snormy(i,j),x(i-1),x(i-1)+abs(u(i,j))*dt,y(j),y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i-1)+abs(u(i,j))*dt,y(j),y(j-1),cnw,csw,cne,cse,vfl(i,j))

			call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse,vol(i,j))

			if(vfl(i,j)<0d0)vfl(i,j)=0d0
			if(vfr(i,j)<0d0)vfr(i,j)=0d0
			if(vol(i,j)<0d0)vol(i,j)=0d0
		  		    
		end if
		  			
	  end if


	end do
  end do
  !solving fluxes for x-y axis

  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i-1,j)/=B.or.mark(i+1,j)/=B) then 

         

	    if(u(i,j)>0d0) then
          
		  vfd=vfr(i,j)
		  vfa=vfl(i+1,j)
		  vd=vol(i,j)
		  va=vol(i+1,j)
		  
		  
        end if
        if(u(i,j)<0d0) then
          
		  vfd=vfl(i+1,j)
		  vfa=vfr(i,j)
		  vd=vol(i+1,j)
		  va=vol(i,j)

        end if

		if(u(i,j)>0d0) then
          sp=slope(i,j)
		  
		end if

        if(u(i,j)<0d0) then
          sp=slope(i+1,j)
         
		end if     
         


        if(abs(sp)>=sc) then   
          
		  vad=va
		  vfad=vfa
          
        else
          
		  vad=vd
		  vfad=vfd
        end if
		
        if(u(i,j)==0d0) then
		  flux_x(i,j)=0d0
		else 
		  cf=max(((abs(u(i,j))*dt*dy-vfad)-(dx*dy-vd)),0d0)
          flux_x(i,j)=min(vfad+cf,vd)         !/dy    
		  !if(flux_x(i,j)<1d-5) flux_x(i,j)=0d0
		  !if(flux_x(i,j)>0.999999) flux_x(i,j)=1d0 
		end if 
	 end if
   end do
 end do
 !solving transport equation for x-sweep
 do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    !fsx(i,j)=(fs(i,j)+1d0/dx*(sign2(u(i-1,j))*flux_x(i-1,j)&
		         ! -sign2(u(i,j))*flux_x(i,j)))/(1d0-dt/dx*(u(i,j)-u(i-1,j)))
        fsx(i,j)=(vol(i,j)+(sign2(u(i-1,j))*flux_x(i-1,j)&
		          -sign2(u(i,j))*flux_x(i,j)))/(1d0-dt/dx*(u(i,j)-u(i-1,j)))    !/(dx*dy)
        fsx(i,j)=fsx(i,j)/(dx*dy)
		!if(fsx(i,j)<0.999d-3) fsx(i,j)=0d0
		!if(fsx(i,j)>0.999) fsx(i,j)=1d0 
	  end if
	end do
  end do




  !y-x sweep
!  write(fptr,'(A)') "y-x sweep"
  slope=0d0
  snormx=0d0
  snormy=0d0
  c=0d0
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

      if(fsx(i,j)<=0d0+1d-5) then	    
		vfl(i,j)=0d0
		vfr(i,j)=0d0
		vol(i,j)=0d0
	  else if(fsx(i,j)>=1d0-1d-5)then
	    vfr(i,j)=abs(v(i,j))*dt*dx !vfr -top
		vfl(i,j)=abs(v(i,j))*dt*dx  !vfl-bottom
		vol(i,j)=dx*dy
	  else
	    if(mark(i,j)/=E.or.mark(i,j)/=F)then
		    call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j)-abs(v(i,j))*dt,cnw,csw,cne,cse,vfr(i,j))

			call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j-1)+abs(v(i,j))*dt,y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j-1)+abs(v(i,j))*dt,y(j-1),cnw,csw,cne,cse,vfl(i,j))

			call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse,vol(i,j))

			if(vfl(i,j)<0d0)vfl(i,j)=0d0
			if(vfr(i,j)<0d0)vfr(i,j)=0d0
			if(vol(i,j)<0d0)vol(i,j)=0d0
		  		    
		end if
		  			
	  end if


	end do
  end do


do i=2,L1-1
 do j=2,M1-1
   vn=vn+vol(i,j)
 end do
end do

do i=2,L1-1
 do j=2,M1-1
   vofsum=vofsum+fsx(i,j)
 end do
end do

!do i=2,L1-1
! do j=2,M1-1
!   fsx(i,j)=fsx(i,j)*(1d0+(vnull-vn)/vofsum*sign2(vnull-vn))
  
!   if(fsx(i,j)<=0d0+1d-5) then	    
!		vfl(i,j)=0d0
!		vfr(i,j)=0d0
!		vol(i,j)=0d0
!	  else if(fsx(i,j)>=1d0-1d-5)then
!	    vfr(i,j)=abs(v(i,j))*dt*dx !vfr -top
!		vfl(i,j)=abs(v(i,j))*dt*dx  !vfl-bottom
!		vol(i,j)=dx*dy
!	  else
!	    if(mark(i,j)/=E.or.mark(i,j)/=F)then
!		    call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse)
!		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j)-abs(v(i,j))*dt,cnw,csw,cne,cse,vfr(i,j))
!
!			call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j-1)+abs(v(i,j))*dt,y(j-1),cnw,csw,cne,cse)
!		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j-1)+abs(v(i,j))*dt,y(j-1),cnw,csw,cne,cse,vfl(i,j))
!
!			call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse)
!		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse,vol(i,j))
!
!			if(vfl(i,j)<0d0)vfl(i,j)=0d0
!			if(vfr(i,j)<0d0)vfr(i,j)=0d0
!			if(vol(i,j)<0d0)vol(i,j)=0d0
!		  		    
!		end if
		  			
!	  end if
! end do
!end do





  


  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i,j-1)/=B.or.mark(i,j+1)/=B) then
	    if(v(i,j)>0d0) then         
		  
		  vfd=vfr(i,j)
		  vfa=vfl(i,j+1)
		  vd=vol(i,j)
		  va=vol(i,j+1)
		  	     
        end if
        if(v(i,j)<0d0) then
          vfd=vfr(i,j+1)
		  vfa=vfl(i,j)
		  vd=vol(i,j+1)
		  va=vol(i,j)     
        end if

		if(v(i,j)>0d0) then
          sp=slope(i,j)
		  
		  
		end if

		if(v(i,j)<0d0) then          
		  sp=slope(i,j+1)
		  
		end if        

		if(abs(sp)>=sc) then   
          vad=va
		  vfad=vfa
		  
        else
          vad=vd
		  vfad=vfd
        end if

        if(v(i,j)==0d0) then
		  flux_y(i,j)=0d0
		else
		  
		   cf=max(((abs(v(i,j))*dt*dx-vfad)-(dx*dy-vd)),0d0)
           flux_y(i,j)=min(vfad+cf,vd)     !/dx
		   !if(flux_y(i,j)<1d-5) flux_y(i,j)=0d0
		   !if(flux_y(i,j)>0.999999) flux_y(i,j)=1d0 
		 end if
	    
	 end if
   end do
 end do
   
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    fs(i,j)=vol(i,j)*(1d0+dt/dy*(v(i,j)-v(i,j-1))) &
		        +(sign2(v(i,j-1))*flux_y(i,j-1)-sign2(v(i,j))*flux_y(i,j)) !/(dx*dy)
        fs(i,j)=fs(i,j)/(dx*dy)
		!if(fs(i,j)<0.999d-3) fs(i,j)=0d0
		!if(fs(i,j)>0.999) fs(i,j)=1d0 
	  end if
   end do
 end do

  slope=0d0
  slope=0d0
  snormx=0d0
  snormy=0d0
  c=0d0
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
else

!y-x sweep
!write(fptr,'(A)') "solving tr.equation y-x sweep,x-y sweep"
  slope=0d0
  snormx=0d0
  snormy=0d0
  c=0d0
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
      

	  if(fs(i,j)<=0d0+1d-5) then	    
		vfl(i,j)=0d0
		vfr(i,j)=0d0
		vol(i,j)=0d0
	  else if(fs(i,j)>=1d0-1d-5)then
	    vfr(i,j)=abs(v(i,j))*dt*dx !vfr -top
		vfl(i,j)=abs(v(i,j))*dt*dx  !vfl-bottom
		vol(i,j)=dx*dy
	  else
	    if(mark(i,j)/=E.or.mark(i,j)/=F)then
		    call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j)-abs(v(i,j))*dt,cnw,csw,cne,cse,vfr(i,j))

			call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j-1)+abs(v(i,j))*dt,y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j-1)+abs(v(i,j))*dt,y(j-1),cnw,csw,cne,cse,vfl(i,j))

			call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse,vol(i,j))

			if(vfl(i,j)<0d0)vfl(i,j)=0d0
			if(vfr(i,j)<0d0)vfr(i,j)=0d0
			if(vol(i,j)<0d0)vol(i,j)=0d0
		  		    
		end if
		  			
	  end if



	end do
  end do

  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i,j-1)/=B.or.mark(i,j+1)/=B) then
	    if(v(i,j)>0d0) then
          vfd=vfr(i,j)
		  vfa=vfl(i,j+1)
		  vd=vol(i,j)
		  va=vol(i,j+1)     
        end if
        if(v(i,j)<0d0) then
          vfd=vfr(i,j+1)
		  vfa=vfl(i,j)
		  vd=vol(i,j+1)
		  va=vol(i,j) 	      
        end if

        if(v(i,j)>0d0) then
          sp=slope(i,j)
          
		end if

		if(v(i,j)<0d0) then
          sp=slope(i,j+1)
		  
		end if

        

        if(abs(sp)>=sc) then   
          vad=va
		  vfad=vfa
        else
          vad=vd
		  vfad=vfd
        end if

        if(v(i,j)==0d0) then
		  flux_y(i,j)=0d0
		else
		  cf=max(((abs(v(i,j))*dt*dx-vfad)-(dx*dy-vd)),0d0)
          flux_y(i,j)=min(vfad+cf,vd)   !/dx
		  !if(flux_y(i,j)<1d-5) flux_y(i,j)=0d0
		  !if(flux_y(i,j)>0.999999) flux_y(i,j)=1d0 
        end if
	 end if
   end do
 end do
   
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    fsy(i,j)=(vol(i,j)+(sign2(v(i,j-1))*flux_y(i,j-1) &
		          -sign2(v(i,j))*flux_y(i,j)))/(1d0-dt/dy*(v(i,j)-v(i,j-1)))  !/(dx*dy)
        fsy(i,j)=fsy(i,j)/(dx*dy)
        !if(fsy(i,j)<0.999d-3) fsy(i,j)=0d0
		!if(fsy(i,j)>0.999) fsy(i,j)=1d0  
	  end if
   end do
 end do

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
      
	  if(fsy(i,j)<=0d0+1d-5) then	    
		vfl(i,j)=0d0
		vfr(i,j)=0d0
		vol(i,j)=0d0
	  else if(fsy(i,j)>=1d0-1d-5)then
	    vfr(i,j)=abs(u(i,j))*dt*dy !vfr -top
		vfl(i,j)=abs(u(i,j))*dt*dy  !vfl-bottom
		vol(i,j)=dx*dy
	  else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!aaaaaa
	    if(mark(i,j)/=E.or.mark(i,j)/=F)then
		    call plc(snormx(i,j),snormy(i,j),x(i)-abs(u(i,j))*dt,x(i),y(j),y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i)-abs(u(i,j))*dt,x(i),y(j),y(j-1),cnw,csw,cne,cse,vfr(i,j))

			call plc(snormx(i,j),snormy(i,j),x(i-1),x(i-1)+abs(u(i,j))*dt,y(j),y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i-1)+abs(u(i,j))*dt,y(j),y(j-1),cnw,csw,cne,cse,vfl(i,j))

			call plc(snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse)
		    call value2(c(i,j),snormx(i,j),snormy(i,j),x(i-1),x(i),y(j),y(j-1),cnw,csw,cne,cse,vol(i,j))

			if(vfl(i,j)<0d0)vfl(i,j)=0d0
			if(vfr(i,j)<0d0)vfr(i,j)=0d0
			if(vol(i,j)<0d0)vol(i,j)=0d0
		  		    
		end if
		  			
	  end if



	end do
  end do 
  !solving fluxes for x-y axis
  do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B.or.mark(i-1,j)/=B.or.mark(i+1,j)/=B) then
	    if(u(i,j)>0d0) then
          vfd=vfr(i,j)
		  vfa=vfl(i+1,j)
		  vd=vol(i,j)
		  va=vol(i+1,j)    
        end if
        if(u(i,j)<0d0) then
          vfd=vfl(i+1,j)
		  vfa=vfr(i,j)
		  vd=vol(i+1,j)
		  va=vol(i,j)	      
        end if

        if(u(i,j)>0d0) then
          sp=slope(i,j)
		  
		end if

        if(u(i,j)<0d0) then
          sp=slope(i+1,j)
          
		end if

        if(abs(sp)>=sc) then   
          vad=va
		  vfad=vfa
        else
		  vad=vd
		  vfad=vfd
        end if
        if(u(i,j)==0d0) then
		  flux_x(i,j)=0d0
		else 
		  cf=max(((abs(u(i,j))*dt*dy-vfad)-(dx*dy-vd)),0d0)
          flux_x(i,j)=min(vfad+cf,vd)      !/dy  
		  !if(flux_x(i,j)<1d-5) flux_x(i,j)=0d0
		  !if(flux_x(i,j)>0.999999) flux_x(i,j)=1d0 
		end if 
	 end if
   end do
 end do
 !solving transport equation for x-sweep
 do i=2,L1-1
    do j=2,M1-1
	  if(mark(i,j)/=B) then
	    fs(i,j)=vol(i,j)*(1d0+dt/dx*(u(i,j)-u(i-1,j)))&
		        +1d0/dx*(sign2(u(i-1,j))*flux_x(i-1,j)-sign2(u(i,j))*flux_x(i,j))
		fs(i,j)=fs(i,j)/(dx*dy)
		!if(fs(i,j)<0.999d-3) fs(i,j)=0d0
		!if(fs(i,j)>0.999) fs(i,j)=1d0 

	  end if
	end do
  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  slope=0d0
  slope=0d0
  snormx=0d0
  snormy=0d0
  c=0d0
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

end if


end subroutine advect_da5