subroutine normal(flag,L1,M1,x,y,ff,mark,snormx,snormy)
use vof2d
implicit none
integer,parameter::n1=2
integer l,m,i,j
logical,intent(in):: flag
integer,intent(in)::L1,M1
real(8),intent(in):: x(ni),y(nj),ff(ni,nj)
integer,intent(in)::mark(ni,nj)
real(8),intent(out):: snormx(ni,nj),snormy(ni,nj)
real(8) a(n1,n1),bb(n1),ss(n1),absnorm,ff2(ni,nj)
real(8) delta, delta1,delta2,xw,xc,xe,ys,yc,yn
real(8) a1(8,2),b1(8)
real(8),parameter::eps=1d-5
real(8),external::sign2
snormx=0d0
snormy=0d0
ff2=0d0


do i=2,L1-1
  do j=2,M1-1
    ff2(i,j)=ff(i,j)
    !print*,ff2(i,j)
  end do
end do
 !выяснить как обрабатывать нормали в полных клетках.    

do i=2,L1-1
  do j=2,M1-1
  !mark(i,j)==S ff(i,j)>0d0.and.ff(i,j)<1d0
    if(ff(i,j)>0d0.and.ff(i,j)<1d0) then
   
    if(mark(i,j)/=B) then

	  yn=y(j+1)+0.5d0*(y(j+2)-y(j+1))
	  yc=y(j)+0.5d0*(y(j+1)-y(j))
	  ys=y(j-1)+0.5d0*(y(j)-y(j-1))

	  xw=x(i-1)+0.5d0*(x(i)-x(i-1))
	  xc=x(i)+0.5d0*(x(i+1)-x(i))
	  xe=x(i+1)+0.5d0*(y(i+2)-y(i+1))

      

       if(flag) then

        if(mark(i+1,j)==B) then
		 l=i+1
		 m=j
		 ff2(i+1,j)=1-ff2(i+i-l,j+j-m)
	    end if

		if(mark(i-1,j)==B) then
		 l=i-1
		 m=j
		 ff2(i-1,j)=1-ff2(i+i-l,j+j-m)
	    end if

		if(mark(i,j+1)==B) then
		 l=i
		 m=j+1
		 ff2(i,j)=1-ff2(i+i-l,j+j-m)
	    end if

		if(mark(i,j-1)==B) then
		 l=i
		 m=j-1
		 ff2(i,j)=1-ff2(i+i-l,j+j-m)
	    end if
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
        if(mark(i+1,j+1)==B) then
		 l=i+1
		 m=j+1
		 ff2(i+1,j+1)=1-ff2(i+i-l,j+j-m)
	    end if
        
		if(mark(i+1,j-1)==B) then
		 l=i+1
		 m=j-1
		 ff2(i+1,j-1)=1-ff2(i+i-l,j+j-m)
	    end if

		if(mark(i-1,j+1)==B) then
		 l=i-1
		 m=j+1
		 ff2(i-1,j+1)=1-ff2(i+i-l,j+j-m)
	    end if
        
		if(mark(i-1,j-1)==B) then
		 l=i-1
		 m=j-1
		 ff2(i-1,j-1)=1-ff2(i+i-l,j+j-m)
	    end if

      end if	

!	  
	  a(1,1)=3d0*(xe-xc)*(xe-xc)+3d0*(xw-xc)*(xw-xc)

	  a(1,2)=(xe-xc)*(yn-yc)+(xw-xc)*(yn-yc)+(xw-xc)*(ys-yc)+(xe-xc)*(ys-yc)	  

	  a(2,1)=(xe-xc)*(yn-yc)+(xw-xc)*(yn-yc)+(xw-xc)*(ys-yc)+(xe-xc)*(ys-yc) !!!!!!!!!

	  a(2,2)=3d0*(yn-yc)*(yn-yc)+3d0*(ys-yc)*(ys-yc)

	  bb(1)=(xe-xc)*(ff2(i+1,j)+ff2(i+1,j+1)+ff2(i+1,j-1)-3d0*ff2(i,j))&
	  +(xw-xc)*(ff2(i-1,j+1)+ff2(i-1,j)+ff2(i-1,j-1)-3d0*ff2(i,j))

	  bb(2)=(yn-yc)*(ff2(i+1,j+1)+ff2(i,j+1)+ff2(i-1,j+1)-3d0*ff2(i,j))&
	  +(ys-yc)*(ff2(i-1,j-1)+ff2(i,j-1)+ff2(i+1,j-1)-3d0*ff2(i,j))
      
      delta=a(1,1)*a(2,2)-a(1,2)*a(2,1)
	  delta1=bb(1)*a(2,2)-bb(2)*a(1,2)
	  delta2=bb(2)*a(1,1)-bb(1)*a(2,1)
	  
	  snormx(i,j)=delta1/delta
      snormy(i,j)=delta2/delta
      

      absnorm=sqrt(snormx(i,j)*snormx(i,j)+snormy(i,j)*snormy(i,j))
	  snormx(i,j)=-snormx(i,j)/absnorm
      snormy(i,j)=-snormy(i,j)/absnorm

      
      if((1d0-abs(snormx(i,j)))<eps) snormx(i,j)=1d0*sign2(snormx(i,j))

      if((1d0-abs(snormy(i,j)))<eps) snormy(i,j)=1d0*sign2(snormy(i,j))

      if(abs(snormx(i,j))<eps) snormx(i,j)=0d0
      
	  if(abs(snormy(i,j))<eps) snormy(i,j)=0d0

!	  print*, snormx(i,j)
!      print*, snormy(i,j)
      ! write(*,*) snormx(i,j)
      !write(*,*) snormy(i,j)
     
	end if
	end if
  end do
end do
end