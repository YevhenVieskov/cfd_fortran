subroutine marking(fb,fs,mark,L1,M1)
use vof2d
implicit none
integer,intent(in)::L1,M1
real(8),intent(in)::fb(ni,nj),fs(ni,nj)
integer,intent(out)::mark(ni,nj)
integer i,j,M2,L2
real(8)::eps=1d-5
L2=L1-1; M2=M1-1
do i=2,L1
 do j=2,M1
  if(fb(i,j)==0d0) mark(i,j)=B
  if(fb(i,j)>0d0.and.fs(i,j)==0d0) mark(i,j)=E
  if(fb(i,j)>0d0.and.fs(i,j)>0d0.and.fs(i,j)<fb(i,j)) mark(i,j)=S
  if(fb(i,j)>0d0.and.fs(i,j)>0d0.and.fs(i,j)==fb(i,j)) mark(i,j)=F
  if(fs(i,j)<=eps)mark(i,j)=E
  if(fs(i,j)>=(1d0-eps))mark(i,j)=F
 end do
end do
end subroutine marking
  
