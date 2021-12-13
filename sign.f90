real(8) function sign2(x)
real(8),intent(in)::x
if(x<0)then
  sign2=-1d0
else
  sign2=1d0
end if
end function sign2