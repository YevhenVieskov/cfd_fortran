c	http://www.pudn.com/ > VOF.rar > superbee.for, change:2001-03-14,size:5391b



c       distance of point(x1,y1) and point(x2,y2) 
      function dist(x1,y1,x2,y2) 
      dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) 
      end  
      function fi(x) 
       fi=max(0.0,min(1.0,2.0*x),min(2.0,x)) 
      end  
      function polate(x1,x2,x3,f1,f2,f3,x)  
      p1=(x-x2)*(x-x3)/((x1-x2)*(x1-x3))*f1 
      p2=(x-x1)*(x-x3)/((x2-x1)*(x2-x3))*f2 
      p3=(x-x1)*(x-x2)/((x3-x1)*(x3-x2))*f3 
      polate=p1+p2+p3 
      end 
c     program main     
c    solve fluid volume function equation using  
c    integrating average Superbee TVD method      
      program main        
      parameter (isize=80,pi=3.14159265)     
      real c0(0:isize,0:isize),c1(0:isize,0:isize) 
      real c2(0:isize,0:isize),gama(0:isize,0:isize) 
      real u(0:isize,0:isize),v(0:isize,0:isize) 
      real x(0:isize),y(0:isize) 
 
c      open(unit=1,file='d:\mytemp\c1.plt') 
c      open(unit=2,file='d:\mytemp\c2.plt') 
c      open(unit=3,file='d:\mytemp\c3.plt') 
c      open(unit=4,file='d:\mytemp\c4.plt') 
       
                   
      write(*,*) 'Input T:' 
      read(*,*) supt 
c      Initialization of parameter dt,dx...       
      h=1.0/isize 
      dt=0.1*h 
      radius=0.4                            
      x0=0.5 
      y0=0.5 
      cki=pi 
c   Initialization of Fluid Volume Function and velocity of fluid       
      do i=0,isize 
         x(i)=h*(i) 
         y(i)=h*(i) 
      enddo          
      do i=0,isize 
        do j=0,isize 
c          u(i,j)=-pi*(y(j)-0.5) 
c          v(i,j)=pi*(x(i)-0.5) 
          u(i,j)=-pi*cos(pi*(x(i)-0.5))*sin(pi*(y(j)-0.5)) 
          v(i,j)=pi*sin(pi*(x(i)-0.5))*cos(pi*(y(j)-0.5))    
c          if(dist(x(i),y(j),x0,y0).lt.radius) then 
c            c0(i,j)=1.0 
c          else 
c            c0(i,j)=0 
c          endif        
c          if(y(j).lt.0.5.and.x(i).gt.0.425.and.x(i).lt.0.575) then 
c            c0(i,j)=0 
c          endif 
c          xx=x(i)-0.5 
c          yy=y(j)-0.5      
c          c0(i,j)=sqrt(xx**2+yy**2)-0.4 
c          if(x(i).ge.0.425.and.x(i).le.0.575.and.y(j).le.0.5) then 
c            if(sqrt(xx**2+yy**2).le.0.4) then 
c              c0(i,j)=min(abs(x(i)-0.425),abs(x(i)-0.575)) 
c            else        
c              temp1=sqrt((x(i)-0.425)**2+(y(j)-0.1)**2) 
c              temp2=sqrt((x(i)-0.575)**2+(y(j)-0.1)**2) 
c              c0(i,j)=min(temp1,temp2) 
c            endif 
c          endif      
          xx=x(i)-0.5 
          yy=y(j)-0.3 
          c0(i,j)=sqrt(xx**2+yy**2)-0.2 
        enddo 
      enddo          
c       Solving... 
      t=0  
      it=0      
      do while(t.lt.supt) 
        if(t+dt.gt.supt) dt=supt-t 
        t=t+dt 
        write(*,*) t                                           
         
        if(t.le.0.5.and.t+dt.gt.0.5)call output(c0,isize,x,y,1) 
        if(t.le.1.0.and.t+dt.gt.1.0)call output(c0,isize,x,y,2)          
        if(t.le.1.5.and.t+dt.gt.1.5)call output(c0,isize,x,y,3) 
        do i=1,isize-1 
          do j=1,isize-1      
            if(abs(c0(i+1,j)-c0(i,j)).gt.1.0e-6) then 
              cita=(c0(i,j)-c0(i-1,j))/(c0(i+1,j)-c0(i,j)) 
            else 
              cita=0 
            endif     
            gama(i,j)=fi(cita)*(c0(i+1,j)-c0(i,j))/h 
          enddo 
        enddo     
        do i=2,isize-2 
          do j=2,isize-2                               
            vx=u(i,j) 
            sigma=dt/h                                            
            d1=-sigma*max(0.,vx)*(c0(i,j)-c0(i-1,j)) 
            d2=sigma/2.0*h*max(0.,vx)*(vx*sigma-1.0) 
     &                              *(gama(i,j)-gama(i-1,j)) 
            d3=-sigma*min(0.,vx)*(c0(i+1,j)-c0(i,j)) 
            d4=sigma/2.0*h*min(0.,vx)*(vx*sigma+1.0) 
     &                              *(gama(i+1,j)-gama(i,j)) 
            c1(i,j)=c0(i,j)+d1+d2+d3+d4 
          enddo 
        enddo 
        do i=1,isize-1 
          do j=1,isize-1  
            if(abs(c1(i,j+1)-c1(i,j)).gt.1.0e-6) then 
               cita=(c1(i,j)-c1(i,j-1))/(c1(i,j+1)-c1(i,j)) 
            else 
               cita=0 
            endif       
            gama(i,j)=fi(cita)*(c1(i,j+1)-c1(i,j))/h 
          enddo 
        enddo     
        do i=2,isize-2 
          do j=2,isize-2                               
            vy=v(i,j) 
            sigma=dt/h                                            
            d1=-sigma*max(0.,vy)*(c1(i,j)-c1(i,j-1)) 
            d2=sigma/2.0*h*max(0.,vy)*(vy*sigma-1.0) 
     &                              *(gama(i,j)-gama(i,j-1)) 
            d3=-sigma*min(0.,vy)*(c1(i,j+1)-c1(i,j)) 
            d4=sigma/2.0*h*min(0.,vy)*(vy*sigma+1.0)         
     &                              *(gama(i,j+1)-gama(i,j)) 
            c2(i,j)=c1(i,j)+d1+d2+d3+d4 
          enddo 
        enddo 
        do i=0,isize 
          do j=0,isize 
            c0(i,j)=c2(i,j) 
          enddo 
        enddo    
        it=it+1  
      enddo 
      call output(c0,isize,x,y,4) 
      end    
        
      subroutine output(c0,isize,x,y,ifile) 
      real c0(0:isize,0:isize),x(0:isize),y(0:isize) 
       
      write(ifile,*)  'TITLE = "EXAMPLE: 3D GEOMETRIES" ' 
      write(ifile,*)  'VARIABLES = "X ", "Y  ", "fi " ' 
      write(ifile,*)  'ZONE T="Floor", I=',isize,' J=',isize,' F=POINT ' 
       
      do 50 i=1,isize 
        do 50 j=1,isize 
          write(ifile,*) x(i),y(j),c0(i,j) 
   50 continue  
      close(ifile) 
      end       

