subroutine vof_vtk(IOpst,filename1,filename2,u,v,p,fs,mark,x,y,nx,ny,c,L1,M1,dx,dy)
use vof2d
implicit none
!integer,parameter:: ni=5,nj=5
integer,intent(in):: IOpst,L1,M1,mark(ni,nj)
real(8),intent(in):: u(ni,nj),v(ni,nj),p(ni,nj),nx(ni,nj),ny(ni,nj),c(ni,nj)
real(8),intent(in)::fs(ni,nj)
real(8),intent(in):: x(ni),y(nj),dx,dy
integer ios
character(200),intent(in)::filename1,filename2 
character(255) fln
integer i,j,ind(ni*nj),in,k,numPoint,k1,L2,M2,L3,M3
real(8) x1(0,0:ni),y1(0,0:nj),xk1,yk1,xk2,yk2,utemp(ni,nj),vtemp(ni,nj)
!x1(0,0:ni)=x(0:ni)
!y1(0,0:nj)=y(0:nj)
!fln=trim(filename)//trim('.vtk')
L2=L1-1

M2=M1-1
L3=L1-2
M3=M1-2
open(IOpst,file=filename1,status='new',iostat=ios) !for pressure,velocity,vof-fraction
!fln=trim(filename)//trim('_surf')//trim('.vtk')
if(free_surface) then
open(IOpst+1,file=filename2,status='new',iostat=ios)
end if

write(IOpst,'(A)')   '# vtk DataFile Version 2.0'
write(IOpst,'(A,A,A)')   'vof result output'
write(IOpst,'(A)')      'ASCII'
!free surface PLIC output
!write(IOpst,'(A)')      'DATASET POLYDATA'
!вычисление количества точек для кусочно-линейного представления свободной поверхности
!вычисление координат точек пересечения кусочно-линейных участков с линиями сетки
!запись в файл количества точек, типа и координат
!определение количества прямых и объема данных
!запись в файл кусочно-линейных отрезков
!запись в файл типа сетки (проверить неструктурированную и неравномерную прямоугольную)
!проинтерполировать скорости из граней ячейки в центры 
!записать скорости как векторы, принадлежащие CELL DATA
!
!записать поля давления  и ф-ции сплошности как CELL DATA, выбрать цветовую карту для сплошности
  
write(IOpst,'(A)')      'DATASET RECTILINEAR_GRID' 
write(IOpst,'(A,1X,I5,1X,I5,1X,I5)')      'DIMENSIONS',L2,M2,1
write(IOpst,'(A,1X,I5,1X,A)')      'X_COORDINATES',L2,'double'  
write(IOpst,'(f8.3)')   X(1:L2)  
write(IOpst,'(A,1X,I5,1X,A)')      'Y_COORDINATES',M2,'double' 
write(IOpst,'(f8.3)')   Y(1:L2)  
write(IOpst,'(A,1X,I5,1X,A)')      'Z_COORDINATES',1,'double'
write(IOpst,'(f8.3)')   0.0 
!point data visualization
write(IOpst,'(A,1X,I5)')   'POINT_DATA', (L2)*(M2)
write(IOpst,'(A,1X,A,1X,A)') 'VECTORS', 'velocity_point', 'double'
utemp=0d0
vtemp=0d0
!utemp(1:L1,1)=u(1:L1,1)
!utemp(1:L1,M1)=u(1:L1,M1)

!vtemp(1,1:M1)=v(1,1:M1)
!vtemp(L1,1:M1)=v(L1,1:M1)

do i=1,L2
  do j=1,M2 !!!
     utemp(i,j)=0.5d0*(u(i,j+1)+u(i,j))
  end do
end do

do i=1,L2 !!!
  do j=1,M2
     vtemp(i,j)=0.5d0*(v(i+1,j)+v(i,j))
  end do
end do



do j=1,M2
  do i=1,L2  
    write(IOpst,'(f8.3,1X,f8.3,1X,f8.3)') utemp(i,j),vtemp(i,j),0.0  
  end do
end do

!cell data visualization
write(IOpst,'(A,1X,I5)')   'CELL_DATA', (L3)*(M3)  


write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'vof_fraction', 'double',1
write(IOpst,'(A,1X,A)') 'LOOKUP_TABLE', 'default'
do j=2,M2
  do i=2,L2  
    write(IOpst,'(f8.3)') fs(i,j)
  end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PLIC surface 
if(free_surface) then
!в оригинальной версии условие fs>0, что приводит к отрицательным значениям индексов, 
!условие fs>1е-5 приводит к положительным индексам, но программа не может построить 
!свободную поверхность 
write(IOpst+1,'(A)')   '# vtk DataFile Version 2.0'
write(IOpst+1,'(A,A,A)')   'PLIC surface'
write(IOpst+1,'(A)')      'ASCII'

write(IOpst+1,'(A)')      'DATASET POLYDATA'
numPoint=0
do j=2,M2
  do i=2,L2 
    if(fs(i,j)>0d0.and.(nx(i,j)/=0d0.or.ny(i,j)/=0d0)) numPoint=numPoint+2
  end do
end do
write(IOpst+1,'(A,1X,I5,1X,A)')      'POINTS',numPoint,'double' 
!solving points of intersection
k=0
do j=2,M2
  do i=2,L2
    if(fs(i,j)>0d0) then

	  if(abs(nx(i,j))==1d0.and.abs(ny(i,j))==0d0) then
	    xk1=(c(i,j)-ny(i,j)*y(j-1))/nx(i,j)
        xk2=(c(i,j)-ny(i,j)*y(j))/nx(i,j)
		k=k+1
		ind(k)=k
		write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') xk1,y(j-1),0.0
		k=k+1
		ind(k)=k
		write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') xk2,y(j),0.0
	  end if

	  if(abs(nx(i,j))==0d0.and.abs(ny(i,j))==1d0) then
	    yk1=(c(i,j)-nx(i,j)*x(i-1))/ny(i,j)
        yk2=(c(i,j)-nx(i,j)*x(i))/ny(i,j)
		k=k+1
		ind(k)=k
		write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') x(i-1),yk1,0.0
		k=k+1
		ind(k)=k
		write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') x(i),yk2,0.0
	  end if
      
	  if(abs(nx(i,j))/=0d0.and.abs(ny(i,j))/=0d0) then
        xk1=(c(i,j)-ny(i,j)*y(j-1))/nx(i,j)
        xk2=(c(i,j)-ny(i,j)*y(j))/nx(i,j)
	    yk1=(c(i,j)-nx(i,j)*x(i-1))/ny(i,j)
        yk2=(c(i,j)-nx(i,j)*x(i))/ny(i,j)

        if(xk1>=x(i-1).and.xk1<=x(i)) then
		  k=k+1
		  ind(k)=k
		  write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') xk1,y(j-1),0.0
		end if
        
		 if(xk2>=x(i-1).and.xk2<=x(i)) then
		  k=k+1
		  ind(k)=k
		  write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') xk2,y(j),0.0
		end if

		if(yk1>=y(j-1).and.yk1<=y(j)) then
		  k=k+1
		  ind(k)=k
		  write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') x(i-1),yk1,0.0
		end if

		if(yk2>=y(j-1).and.yk2<=y(j)) then
		  k=k+1
		  ind(k)=k
		  write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') x(i),yk2,0.0
		end if		
	  end if
	end if
  end do 
end do

k=0
do j=1,M1
  do i=1,L1  
    if(fs(i,j)>0d0.and.(nx(i,j)/=0d0.or.ny(i,j)/=0d0)) k=k+1
  end do
end do

write(IOpst+1,'(A,1X,I5,1X,I5)')      'LINES',k,k*3 

do k=1,numPoint,2
  !if(k==1)then
  !  k1=k
  !else 
   ! k1=k+1
  !end if
  write(IOpst+1,'(I5,1X,I5,1X,I5)') 2,ind(k)-1,ind(k+1)-1
end do

end if

end subroutine vof_vtk




subroutine vof_vtk2(IOpst,filename1,filename2,u,v,p,fs,mark,x,y,nx,ny,c,L1,M1,dx,dy)
use vof2d
implicit none
!integer,parameter:: ni=5,nj=5
integer,intent(in):: IOpst,L1,M1,mark(ni,nj)
real(8),intent(in):: u(ni,nj),v(ni,nj),p(ni,nj),nx(ni,nj),ny(ni,nj),c(ni,nj)
real(8),intent(in)::fs(ni,nj)
real(8),intent(in):: x(ni),y(nj),dx,dy
integer ios
character(200),intent(in)::filename1,filename2 
character(255) fln
integer i,j,ind(ni*nj),in,k,numPoint,k1
real(8) x1(1,1:ni),y1(1,1:nj),xk1,yk1,xk2,yk2,utemp(ni,nj),vtemp(ni,nj)
x1(1,1:ni)=x(1:ni)
y1(1,1:nj)=y(1:nj)
!fln=trim(filename)//trim('.vtk')
open(IOpst,file=filename1,status='new',iostat=ios) !for pressure,velocity,vof-fraction
!fln=trim(filename)//trim('_surf')//trim('.vtk')
if(free_surface) then
open(IOpst+1,file=filename2,status='new',iostat=ios)
end if

write(IOpst,'(A)')   '# vtk DataFile Version 2.0'
write(IOpst,'(A,A,A)')   'vof result output'
write(IOpst,'(A)')      'ASCII'
!free surface PLIC output
!write(IOpst,'(A)')      'DATASET POLYDATA'
!вычисление количества точек для кусочно-линейного представления свободной поверхности
!вычисление координат точек пересечения кусочно-линейных участков с линиями сетки
!запись в файл количества точек, типа и координат
!определение количества прямых и объема данных
!запись в файл кусочно-линейных отрезков
!запись в файл типа сетки (проверить неструктурированную и неравномерную прямоугольную)
!проинтерполировать скорости из граней ячейки в центры 
!записать скорости как векторы, принадлежащие CELL DATA
!
!записать поля давления  и ф-ции сплошности как CELL DATA, выбрать цветовую карту для сплошности
  
write(IOpst,'(A)')      'DATASET RECTILINEAR_GRID' 
write(IOpst,'(A,1X,I5,1X,I5,1X,I5)')      'DIMENSIONS',L1,M1,1
write(IOpst,'(A,1X,I5,1X,A)')      'X_COORDINATES',L1,'double'  
write(IOpst,'(f8.3)')   X1(1,1:L1)  
write(IOpst,'(A,1X,I5,1X,A)')      'Y_COORDINATES',M1,'double' 
write(IOpst,'(f8.3)')   y1(1,1:M1)  
write(IOpst,'(A,1X,I5,1X,A)')      'Z_COORDINATES',1,'double'
write(IOpst,'(f8.3)')   0.0 
!point data visualization
write(IOpst,'(A,1X,I5)')   'POINT_DATA', (L1)*(M1)
write(IOpst,'(A,1X,A,1X,A)') 'VECTORS', 'velocity_point', 'double'

utemp(1,1)=u(1,1); vtemp(1,1)=v(1,1) !left bottom
utemp(1,M1)=u(1,M1); vtemp(1,M1)=v(1,M1) !left top

utemp(L1,1)=u(L1,1); vtemp(L1,1)=0d0 !right bottom
utemp(L1,M1)=u(L1,M1); vtemp(L1,M1)=v(L1,M1) !right top

do j=2,M1-1
  utemp(1,j)=0.5d0*(u(1,j)+u(1,j+1))!left boundary
  vtemp(1,j)=0d0
end do

do j=2,M1-1
  utemp(L1,j)=0.5d0*(u(L1,j)+u(L1,j+1)) !right boundary
  vtemp(L1,j)=0d0
end do

do i=2,L1-1
  utemp(i,1)=0d0 !bottom boundary
  vtemp(i,1)=0.5d0*(v(i,1)+v(i+1,1))
end do

do i=2,L1-1
  utemp(i,M1)=0d0 !top boundary
  vtemp(i,M1)=0.5d0*(v(i,M1)+v(i+1,M1))
end do

do j=2,M1-1
  do i=2,L1-1  
     utemp(i,j)=0.5d0*(u(i,j)+u(i+1,j))
	 vtemp(i,j)=0.5d0*(v(i,j)+v(i+1,j))
  end do
end do

do j=1,M1
  do i=1,L1  
    write(IOpst,'(f8.3,1X,f8.3,1X,f8.3)') utemp(i,j),vtemp(i,j),0.0  
  end do
end do
!write(IOpst,'(f8.3,1X,f8.3,1X,f8.3)') 0.5d0*(u(i,j)+u(i+1,j)),0.5d0*(v(i,j)+v(i+1,j)),0.0
!cell data visualization
write(IOpst,'(A,1X,I5)')   'CELL_DATA', (L1-1)*(M1-1)  
write(IOpst,'(A,1X,A,1X,A)') 'VECTORS', 'velocity', 'double'
do j=2,M1
  do i=2,L1  
    write(IOpst,'(f8.3,1X,f8.3,1X,f8.3)') 0.5d0*(u(i-1,j)+u(i,j)),0.5d0*(v(i,j-1)+v(i,j)),0.0
  end do
end do

write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'velocity_magnitude', 'double',1
write(IOpst,'(A,1X,A)') 'LOOKUP_TABLE', 'default'
do j=2,M1
  do i=2,L1  
    write(IOpst,'(f8.3)') sqrt(u(i,j)**2+v(i,j)**2)
  end do
end do


write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'pressure', 'double',1
write(IOpst,'(A,1X,A)') 'LOOKUP_TABLE', 'default'
do j=2,M1  
  do i=2,L1  
    write(IOpst,'(f20.3)') p(i,j)
  end do
end do

write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'vof_fraction', 'double',1
write(IOpst,'(A,1X,A)') 'LOOKUP_TABLE', 'default'
do j=2,M1
  do i=2,L1  
    write(IOpst,'(f8.3)') fs(i,j)
  end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PLIC surface 
if(free_surface) then
write(IOpst+1,'(A)')   '# vtk DataFile Version 2.0'
write(IOpst+1,'(A,A,A)')   'PLIC surface'
write(IOpst+1,'(A)')      'ASCII'

write(IOpst+1,'(A)')      'DATASET POLYDATA'
numPoint=0
do j=2,M1
  do i=2,L1 
    if(mark(i,j)==S) numPoint=numPoint+2
  end do
end do
write(IOpst+1,'(A,1X,I5,1X,A)')      'POINTS',numPoint,'double' 
!solving points of intersection
k=0
do j=2,M1
  do i=2,L1
    if(mark(i,j)==S) then

	  if(abs(nx(i,j))==1d0.and.abs(ny(i,j))==0d0) then
	    xk1=(c(i,j)-ny(i,j)*y(j-1))/nx(i,j)
        xk2=(c(i,j)-ny(i,j)*y(j))/nx(i,j)
		k=k+1
		ind(k)=k
		write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') xk1,y(j-1),0.0
		k=k+1
		ind(k)=k
		write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') xk2,y(j),0.0
	  end if

	  if(abs(nx(i,j))==0d0.and.abs(ny(i,j))==1d0) then
	    yk1=(c(i,j)-nx(i,j)*x(i-1))/ny(i,j)
        yk2=(c(i,j)-nx(i,j)*x(i))/ny(i,j)
		k=k+1
		ind(k)=k
		write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') x(i-1),yk1,0.0
		k=k+1
		ind(k)=k
		write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') x(i),yk2,0.0
	  end if
      
	  if(abs(nx(i,j))/=0d0.and.abs(ny(i,j))/=0d0) then
        xk1=(c(i,j)-ny(i,j)*y(j-1))/nx(i,j)
        xk2=(c(i,j)-ny(i,j)*y(j))/nx(i,j)
	    yk1=(c(i,j)-nx(i,j)*x(i-1))/ny(i,j)
        yk2=(c(i,j)-nx(i,j)*x(i))/ny(i,j)

        if(xk1>=x(i-1).and.xk1<=x(i)) then
		  k=k+1
		  ind(k)=k
		  write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') xk1,y(j-1),0.0
		end if
        
		 if(xk2>=x(i-1).and.xk2<=x(i)) then
		  k=k+1
		  ind(k)=k
		  write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') xk2,y(j),0.0
		end if

		if(yk1>=y(j-1).and.yk1<=y(j)) then
		  k=k+1
		  ind(k)=k
		  write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') x(i-1),yk1,0.0
		end if

		if(yk2>=y(j-1).and.yk2<=y(j)) then
		  k=k+1
		  ind(k)=k
		  write(IOpst+1,'(f8.3,1X,f8.3,1X,f8.3)') x(i),yk2,0.0
		end if		
	  end if
	end if
  end do 
end do

k=0
do j=2,M1
  do i=2,L1  
    if(mark(i,j)==S) k=k+1
  end do
end do

write(IOpst+1,'(A,1X,I5,1X,I5)')      'LINES',k,k*3 

do k=1,numPoint,2
  !if(k==1)then
  !  k1=k
  !else 
   ! k1=k+1
  !end if
  write(IOpst+1,'(I5,1X,I5,1X,I5)') 2,ind(k)-1,ind(k+1)-1
end do

end if

end subroutine vof_vtk2

