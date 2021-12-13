 subroutine output(IOpst,u,v,p,fs,mark,x,y,nx,ny,c,L1,M1,dx,dy,num)      

use vof2d
implicit none       
integer,intent(in):: IOpst,L1,M1,mark(ni,nj),num
real(8),intent(in):: u(ni,nj),v(ni,nj),p(ni,nj),nx(ni,nj),ny(ni,nj),c(ni,nj)
real(8),intent(in)::fs(ni,nj)
real(8),intent(in):: x(ni),y(nj),dx,dy
integer ios

      
	character(50)	:: szNumber
      character(41)	:: szFileName
	character(41)	:: szFileNameV
	character(41)	:: szFileNameP
	character(41)	:: szFileNameVOF
    character(200) filename1,filename2
	integer fptr,i,j,is,ie,M2,L2
!      generation rezult files v zadannij moment vremeni
      if(num .ne. 0) then
	 Write(szNumber,'(i6.6)') num
	!szFileName= Trim("../voffunc2/data/")//Trim("uvf")//Trim(szNumber)//Trim(".dat")

    !szFileNameV= Trim("../voffunc/data/")//Trim("v")//Trim(szNumber)//Trim(".dat")

    !szFileNameP= Trim("../voffunc/data/")//Trim("p")//Trim(szNumber)//Trim(".dat")

	!szFileNameVOF= Trim("../voffunc/data/")//Trim("vof")//Trim(szNumber)//Trim(".dat")


L2=L1-1;M2=M1-1
      fptr=IOpst
     ! open( fptr,file=szFileName,status='new')	  
      !call PRNT(u,"u",fptr,L1,M1)
	  !write(fptr,*) "   "
	  !call PRNT(v,"v",fptr,L1,M1)
	  !write(fptr,*) "   "
	  !call PRNT(fs,"vof",fptr,L1,M1)
	  !write(fptr,*) "   "
        


	  fptr=fptr+1
	  filename1=Trim("../voffunc3/vtk/")//trim('uvpf')//Trim(szNumber)//trim('.vtk')



     filename2=Trim("../voffunc3/vtk/")//trim('surf')//Trim(szNumber)//trim('.vtk')

    call vof_vtk(fptr,filename1,filename2,u,v,p,fs,mark,x,y,nx,ny,c,L1,M1,dx,dy)


	end if
      

end subroutine output


!SUBROUTINE PRNT(FI,TITLE,numfile,NNI,NNJ)
!use vof2d
!real(8),intent(in)::FI(ni,nj)    
!CHARACTER*6,intent(in):: TITLE
!integer,intent(in)::NNI,NNJ
!      WRITE(numfile,20) TITLE 
!      IS=-11
!  100 IS=IS+12
!      IE=IS+11
!      IE=MIN(NNI,IE)
!      WRITE(numfile,21) (I,I=IS,IE) 
!      WRITE(numfile,22) 
!      DO J=NNJ,1,-1
!        WRITE(numfile,23) J,(FI(I,J),I=IS,IE) 
!      END DO
!      IF(IE.LT.NNI) GO TO 100
!   20 FORMAT(2X,26('*-'),7X,A6,7X,26('-*')) 
!   21 FORMAT(3X,'I = ',I3,11I10)
!   22 FORMAT(2X,'J')
!   23 FORMAT(1X,I3,1P12E10.2) 
!      RETURN
!      END  