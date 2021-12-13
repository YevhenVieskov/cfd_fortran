subroutine voftvd
!fi=1,ak=0.5
do i=2,im1
 do j=2,jm1
    dfl(i,j)=minmod((fs(i,j)-fs(i-1,j)),4d0*(fs(i+1,j)-fs(i,j)))
    dfr(i,j)=minmod((fs(i+1,j)-fs(i,j)),4d0*(fs(i,j)-fs(i-1,j)))
end do
end do 

do i=2,im1
do j=2,jm1
 
end do
end do

do i=2,im1
 do j=2,jm1   
if(u(i,j)<=0d0) then
  flux_x(i,j)=fr(i,j)*u(i,j)
else
  flux_x(i,j)=fl(i,j)*u(i,j)
end if
end do
end do
end

real(8) function minmod(a,b)
real(8) a,b
if(a*b<=0d0)then
  minmod=0d0
else
 minmod=sign(a,1)*min(abs(a),abs(b))
end if
end