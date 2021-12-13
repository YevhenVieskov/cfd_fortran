subroutine plc(nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse)
real(8),intent(in):: nx,ny,xw,xe,yn,ys
real(8),intent(out):: cnw,csw,cne,cse

!real(8) cnw1,csw1,cne1,cse1
cnw=nx*xw+ny*yn
csw=nx*xw+ny*ys
cne=nx*xe+ny*yn
cse=nx*xe+ny*ys

end

subroutine value(c,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,vof)
real(8) nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,vof,c
real(8) ai1,ai2,ai3,ai4,v,vf,ai22,ai23,ai42,ai43,x,y 

if(nx==1d0.and.ny==0d0) then
  v=(yn-ys)*(c/nx-xw)


else if(nx==-1d0.and.ny==0d0) then
  v=(yn-ys)*(xe-c/nx) 

else if(nx==0d0.and.ny==1d0) then
 v=(xe-xw)*(c/ny-ys)

else if(nx==0d0.and.ny==-1d0) then
 v=(xe-xw)*(yn-c/ny)
else

 ai1=abs((min(cnw,c)-min(csw,c))/(max(cnw,csw)-min(cnw,csw)))*(yn-ys)

 ai2=abs((min(cse,c)-min(csw,c))/(max(cse,csw)-min(cse,csw)))*(xe-xw)

 ai3=abs((min(cne,c)-min(cse,c))/(max(cne,cse)-min(cne,cse)))*(yn-ys)

 ai4=abs((min(cne,c)-min(cnw,c))/(max(cne,cnw)-min(cne,cnw)))*(xe-xw)

 v=-(xw-c*nx)/2d0*ai1-(ys-c*ny)/2d0*ai2+(xe-c*nx)/2d0*ai3+(yn-c*ny)/2d0*ai4

end if

vf=(yn-ys)*(xe-xw)
vof=v/vf
end

subroutine value2(c,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,v)
real(8) nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,c
real(8) ai1,ai2,ai3,ai4,v,vf,ai22,ai23,ai42,ai43,x,y 

if(nx==1d0.and.ny==0d0) then
  v=(yn-ys)*(c/nx-xw)


else if(nx==-1d0.and.ny==0d0) then
  v=(yn-ys)*(xe-c/nx) 

else if(nx==0d0.and.ny==1d0) then
 v=(xe-xw)*(c/ny-ys)

else if(nx==0d0.and.ny==-1d0) then
 v=(xe-xw)*(yn-c/ny)
else

 ai1=abs((min(cnw,c)-min(csw,c))/(max(cnw,csw)-min(cnw,csw)))*(yn-ys)

 ai2=abs((min(cse,c)-min(csw,c))/(max(cse,csw)-min(cse,csw)))*(xe-xw)

 ai3=abs((min(cne,c)-min(cse,c))/(max(cne,cse)-min(cne,cse)))*(yn-ys)

 ai4=abs((min(cne,c)-min(cnw,c))/(max(cne,cnw)-min(cne,cnw)))*(xe-xw)

 v=-(xw-c*nx)/2d0*ai1-(ys-c*ny)/2d0*ai2+(xe-c*nx)/2d0*ai3+(yn-c*ny)/2d0*ai4

end if


end


subroutine reset_xy(xw,xe,yn,ys,fs1,fb1,aw,ae,an,as,xe_new,xw_new)
logical a1,a2,a3,a4
real(8) xw,xe,yn,ys,fs1,fb1,aw,ae,an,as,xl,vcell,vsolid

a1=ae<1d0.and.an<1d0.and.((aw==1d0.and.as==1d0).or.(aw==0d0.and.as==0d0))

a2=as<1d0.and.ae<1d0.and.((aw==1d0.and.an==1d0).or.(aw==0d0.and.an==0d0))

a3=as<1d0.and.aw<1d0.and.((as==1d0.and.ae==1d0).or.(as==0d0.and.ae==0d0))

a4=an<1d0.and.aw<1d0.and.((as==1d0.and.ae==1d0).or.(as==0d0.and.ae==0d0))

vcell=(yn-ys)*(xe-xw)
vsolid=(1d0-fb1)*vcell
xl=vsolid/(yn-ys)
if(a1.or.a2) then
  xe_new=xw+xl
  xw_new=xw
end if
if(a3.or.a4) then
  xw_new=xe-xl
  xe_new=xe
end if
end


subroutine planec(fs1,nx,ny,xw,xe,yn,ys,c)
real(8),intent(in)::fs1,nx,ny,xw,xe,yn,ys
real(8),intent(out)::c
real(8) eps,cnw,csw,cne,cse,vofc,cm
real(8) vofnw,vofsw,vofne,vofse
eps=1d-5
call plc(nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse)
call value(cnw,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,vofnw)
call value(csw,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,vofsw)
call value(cne,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,vofne)
call value(cse,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,vofse)

if(fs1==vofsw) then
  c=csw
  return
end if

if(fs1==vofne) then
  c=cne
  return
end if

if(fs1==vofse) then
  c=cse
  return
end if

if(fs1==vofnw) then
  c=cnw
  return
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(nx<0d0.and.ny>0d0)then
 if(vofsw<vofne) then
  if(fs1>vofse.and.fs1<vofsw) then
     cl=cse
     cu=csw
   end if

   if(fs1>vofsw.and.fs1<vofne) then
      cl=csw
      cu=cne
    end if

    if(fs1>vofne.and.fs1<vofnw) then
      cl=cne
      cu=cnw
    end if
 end if

 if(vofsw>vofne) then
  if(fs1>vofse.and.fs1<vofne) then
     cl=cse
     cu=cne
   end if
   
   if(fs1>vofne.and.fs1<vofsw) then
      cl=cne
      cu=csw
    end if

    if(fs1>vofsw.and.fs1<vofnw) then
      cl=csw
      cu=cnw
    end if
 end if

 if(vofsw==vofne) then
   if(fs1<vofne) then
     cl=cse
	 cu=cne
   else if(fs1>vofne)then
     cl=cne
	 cu=cnw     
   end if
 end if

end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(nx>0d0.and.ny<0d0)then
 if(vofne<vofsw) then
  if(fs1>vofnw.and.fs1<vofne) then
     cl=cnw
     cu=cne
   end if

   if(fs1>vofne.and.fs1<vofsw) then
      cl=cne
      cu=csw
    end if

    if(fs1>vofsw.and.fs1<vofse) then
      cl=csw
      cu=cse
    end if
 end if

 if(vofne>vofsw) then
  if(fs1>vofnw.and.fs1<vofsw) then
     cl=cnw
     cu=csw
   end if

   if(fs1>vofsw.and.fs1<vofne) then
      cl=csw
      cu=cne
    end if

    if(fs1>vofne.and.fs1<vofse) then
      cl=csw
      cu=cse
    end if
 end if

 if(vofne==vofsw) then
   if(fs1<vofne) then
     cl=cnw
	 cu=cne
   else if(fs1>vofne)then
     cl=cne
	 cu=cse     
   end if
 end if


end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(nx>0d0.and.ny>0d0)then   !!!!!!!
 if(vofse<vofnw) then
  if(fs1>vofsw.and.fs1<vofse) then
     cl=csw
     cu=cse
   end if

   if(fs1>vofse.and.fs1<vofnw) then
      cl=cse
      cu=cnw
    end if

    if(fs1>vofnw.and.fs1<vofne) then
      cl=cnw
      cu=cne
    end if
  end if

if(vofse>vofnw) then !!
  if(fs1>vofsw.and.fs1<vofnw) then
     cl=csw
     cu=cnw
   end if

   if(fs1>vofnw.and.fs1<vofse) then
      cl=cnw
      cu=cse
    end if

    if(fs1>vofse.and.fs1<vofne) then
      cl=cse
      cu=cne
    end if
  end if
  
  if(vofnw==vofse) then
   if(fs1<vofse) then
     cl=csw
	 cu=cse
   else if(fs1>vofse)then
     cl=cse
	 cu=cne     
   end if
 end if


end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(nx<0d0.and.ny<0d0)then
 if(vofnw<vofse) then
  if(fs1>vofne.and.fs1<vofnw) then
     cl=cne
     cu=cnw
   end if

   if(fs1>vofnw.and.fs1<vofse) then
      cl=cnw
      cu=cse
    end if

    if(fs1>vofse.and.fs1<vofsw) then
      cl=cse
      cu=csw
    end if
  end if

  if(vofnw>vofse) then
  if(fs1>vofne.and.fs1<vofse) then
     cl=cne
     cu=cse
   end if

   if(fs1>vofse.and.fs1<vofnw) then
      cl=cse
      cu=cnw
    end if

    if(fs1>vofnw.and.fs1<vofsw) then
      cl=cnw
      cu=csw
    end if
  end if
  
  if(vofnw==vofse) then
   if(fs1<vofse) then
     cl=cne
	 cu=cse
   else if(fs1>vofse)then
     cl=cse
	 cu=csw     
   end if
 end if

end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(nx==0d0.and.ny==1d0)then
  if(fs1>vofsw.and.fs1<vofnw) then
     cl=csw
     cu=cnw
   end if
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(nx==0d0.and.ny==-1d0)then
  if(fs1>vofnw.and.fs1<vofsw) then
     cl=cnw
     cu=csw
   end if
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(nx==1d0.and.ny==0d0)then
  if(fs1>vofsw.and.fs1<vofse) then
     cl=csw
     cu=cse
   end if
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(nx==-1d0.and.ny==0d0)then
  if(fs1>vofse.and.fs1<vofsw) then
     cl=cse
     cu=csw
   end if
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(nx==0d0.and.ny==0d0)then
  c=0d0
  return
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vofc=0d0

do while(abs(vofc-fs1)>eps)
  cm=0.5d0*(cl+cu)
  call value(cm,nx,ny,xw,xe,yn,ys,cnw,csw,cne,cse,vofc)
  if(vofc<fs1) cl=cm
  if(vofc>fs1) cu=cm
  if(vofc==fs1) c=cm
end do

c=cm
return
end

