module vof2d
integer,parameter::ni=200,nj=200
integer,parameter::B=0,E=1,S=2,F=3
integer,parameter::nd=4,na=2
integer,parameter:: id=ni*nd,jd=nj*nd
integer,parameter:: ia=ni*na,ja=nj*na
!integer M1,M2,L1,L2
!integer ind
real(8):: pi=3.14159,g=9.81,big=1d30,small=1d-30
!глобальные настройки расчета
logical dfd,const_ex, mass_cons,stension,pinterp,da,yangs,lhf,free_surface,cd,ud
real(8) vnull
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!real(8) xmax,xmin,ymax,ymin,XL,YL
!real(8) hx,hy
!real(8) x(ni),y(nj),ax(ni,nj),ay(ni,nj),fb(ni,nj),fs(ni,nj)
!real(8) xcc(ni,nj),ycc(ni,nj),u(ni,nj),v(ni,nj),p(ni,nj)
!real(8) xsub(nd*ni),ysub(nd*nj),hxsub,hysub
!real(8) h(ni,nj),curv(ni,nj) !,mark(ni,nj)
!real(8) rho,amu,anu,theta,gx,gy
!real(8) snormx(ni,nj),snormy(ni,nj)
!properties
!logical dfd, const_ex, lin_ex, mass_cons
!dfd-определ€ет тип аппроксимации диффуных членов,const_ex-константна€ экстрапол€ци€
!lin_ex-линейна€ экстрапол€ци€, mass_cons-сохранение массы  
!ni,nj-максимальные размеры массивов
!L1,M1-верхние границы сеточных линий в направлении x и y соответственно
!E-пуста€ €чейка
!B-€чейка, полностью наход€ща€с€ в твердом теле
!F-€чейка, заполненна€ жидкостью
!S-€чейка, содержаща€ свободную поверхность
!XL,YL-размеры пр€моугольной сетки 
!hx,hy-размеры €чейки
!L2,M2-
!x,y-координаты узлов сетки
!xc,yc-координаты центров €чеек сетки
!ax,ay-апертуры граней €чейки
!fb-объемные апертуры
!fs-апертуры поверхности
!h-значени€ локальных функций высоты
!curv-значени€ кривизны
!mark-маркеры €чеек
!nd-параметр,задающий число узлов подсетки
!xsub,ysub-координаты центров €чеек подсетки
!hxsub,hysub-размеры €чеек подсетки
contains
subroutine setup()
 dfd=.false.
 const_ex=.false.
 mass_cons=.true.
 stension=.false.
 s_angle=.false.
 d_angle=.false.
 pinterp=.false.
 da=.true.
 yangs=.false.
 free_surface=.true.
 cd=.false.
 ud=.true.
end subroutine setup
end module vof2d