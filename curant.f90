subroutine curant(u,v,L1,M1,dx,dy,dt,cfl,t)
use vof2d
implicit none
integer,intent(in)::L1,M1
real(8),intent(in)::t,dx,dy,cfl,u(ni,nj),v(ni,nj)
real(8),intent(out)::dt
integer i,j
real(8) umax,vmax,du,dv
if (cfl>=1d-10) then
  umax=1d-10
  vmax=1d-10
  do i=1,L1
    do j=1,M1
	  if(abs(u(i,j))>umax) umax=abs(u(i,j))
	end do
  end do

  do i=1,L1
    do j=1,M1
      if(abs(v(i,j))>vmax) vmax=abs(v(i,j))
	end do
  end do

  du=dx/umax
  dv=dy/vmax
  if(du<dv)then
   dt=du
  else
   dt=dv
  end if
   !dt=0.1d0
  dt=cfl*dt
end if
end subroutine curant


