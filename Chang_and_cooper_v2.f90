!call system_clock(t) !计算时间
!本程序计算电子谱的演化2013.10.20
!用Chang & Cooper1970的方法
!对第一版做了优化
!2014.10.16加入了自恰的解出Uph和gamma_c


Module constant
implicit none
real*8, parameter:: pi=3.14159265359
real*8, parameter:: c_speed=2.9979d10           !光速
real*8, parameter:: q_e=4.8032d-10         !电子电量
real*8, parameter:: m_e=9.11d-28           !电子质量
real*8, parameter:: sigma_T=6.65d-25      !汤姆森截面
real*8, parameter:: h_planck=6.626d-27           !普朗克常数
real*8, parameter:: re=2.818d-13          !电子半径
real*8, parameter:: m_p=1.673d-24          !质子质量 
end module


        Module Ng       !用于Nga函数的插值
        implicit none
        integer,parameter::n=100
        real*8 x(n-1),uu(n-1)
        end module

Module initial_parameter !包括磁场，R等
implicit none
real*8 B_m,Radius,p,gamma_low,gamma_up_inj,gamma_low_inj,gamma_bulk
end module


	module comptonY
		real*8 comp_Y
	end module

      Module gaus_int
	use initial_parameter
      implicit none
      integer,parameter::nnn=39               !13*5可以调整，越大越精确
!	integer nnn
	real*8 ww1(nnn),ww2(nnn),ww3(nnn),xxx(nnn),nuuu(nnn),ttt(nnn),nu1,nu2      
      end module

program main
use Ng
use initial_parameter
use gaus_int
use comptonY
use constant
implicit none
!integer,parameter::n=1000

real*8 aa(n-1),cc(n-1),bb(n-1),C1(n-1),C2(n-1),C3(n-1),A_A(n-1),r(n-1)      !,x(n),uu(n)
real*8 A,B,C,delta,Q,u   !外部函数
real*8 dx,dt,ratio,deltax,dd,t_tot,ab,inj_t,t_va,Inuu,Inu,rr,L_gamma
real*8 Bi05,Ci05,Bi_05,Ci_05,w05,w_05,delta05,delta_05,loggamma_low,time,r0,q0,B_m0,E_tot,qa0,time1
integer i,j,kk,kabs
character(len=9)::Bm,tot,YY,lumi,Radii,qq0
character(len=55)::nametmp1,nametmp
external a,b,c,delta,u,Q,Inu
common /dx/ dx       !w函数里需要用到dx

common /j/ j        !初始时刻
common /tim/ time,r0
common /dtt/ dt
common /qq0/  q0
common /abb/ ab
common /rr0/ rr
common /kabss/ kabs
common /times/ t_tot
!open(3,file='Changcooper_v2_Y1_1.txt')



p=2.9
gamma_low=3.d0                   !这里是电子的最低能量，但不是注入的最低能量
gamma_low_inj=1.d3                !注入最低能量 
gamma_up_inj=1.d5		  !注入最高能量
nu1=1.d9                          !光子谱上下限
nu2=1.d19                         !
comp_Y=0.0                         !Compton因子

!kq0=1.d15                              !改正因子
gamma_bulk=300.
ab=0.d0                           !=1考虑绝热,=0不考虑
t_va=0.1
inj_t=t_va*gamma_bulk            !t_va观测到的光变时标，inj_t为共动系注入时间,一般可以取为共动系动力学时间

r0=1.d16      !initial radius

Radius=r0/gamma_bulk   !这里应该是共动系的壳层尺度,如果取为0.，则没有自吸收


kabs=1   !k=0不考虑自吸收



t_tot=2.d1                       !计算到的时间
B_m0=1.d5
L_gamma=1.d54
!dt=1.d-2
dt=t_tot/10.
!q0=L_gamma/m_p/c_speed**2*(p-1.)/gamma_low_inj**(1.-p)/(4.*pi*r0**3) !1.d63/1.d15!电子谱归一化系数

q0=1.d30
!print*,q0/qa0
!print*,q0
!pause
!内激波内能为n*mp*c^2
deltax=100.
write(Bm,'(d8.1)')B_m0
write(tot,'(d9.2)')t_tot
write(YY,'(f4.1)')comp_Y
!write(Lumi,'(d9.2)')L_gamma
write(qq0,'(d9.2)')q0
write(Radii,'(d9.2)')r0
nametmp='q0_'//trim(adjustl(adjustr(qq0)))//'Y_'//trim(adjustl(adjustr(YY)))//'B_'//trim(adjustl(adjustr(Bm)))//'t_'//&
&trim(adjustl(adjustr(tot)))//'r0_'//trim(adjustl(adjustr(Radii)))//'.txt'
!nametmp1='f'//trim(adjustl(adjustr(lumi)))//trim(adjustl(adjustr(YY)))//trim(adjustl(adjustr(Bm)))//trim(adjustl(adjustr(tot)))//'.txt'

print*, trim(adjustl(adjustr(nametmp)))
!print*,'g_r_7c'//trim(adjustl(adjustr(nametmp)))
open(3,file='g1_'//trim(adjustl(adjustr(nametmp))))
open(4,file='f1_'//trim(adjustl(adjustr(nametmp))))

!print*,(3,file='sdaaaaaaaaaaaaaaaaaaaaaaaaaaddddddddddddddddddddddddddddddddddd.txt')


aa(1)=0.d0
cc(n-1)=0.d0

deltax=dlog(gamma_up_inj)-dlog(gamma_low)

dx=deltax/n             !(ln(10**4)-ln3)/n,从洛伦兹因子3（对数x=1.多）算到1万（对数9.2）
loggamma_low=dlog(gamma_low)

kk=1            !kk=0 for delta function like injection and 1 for contineous one. 


ratio=dx/dt

!dd=deltax/(1.d0*n) !deltax是要计算的x的总长度，如从x=0计算到x=100,则deltax=100
!print*,nnn
!nnn=65
call gaus_point(gamma_low,gamma_up_inj,xxx,ww1,nnn)   !log10 Lorentz factors of electrons, used for the following integration                                                  

call gaus_point(nu1,nu2,nuuu,ww2,nnn)  !frequency integration points
  
call gaus_point(0.d0,t_tot,ttt,ww3,nnn) !integration over time 

loop2: do j=1,int(t_tot/dt)      !10是想要算到的时间 

time=dt*j
time1=time
rr=c_speed*time*gamma_bulk+r0
!print*,time,rr
if (ab==1.d0) then
	B_m=B_m0*((time+r0/(c_speed*gamma_bulk))/(r0/(c_speed*gamma_bulk)))**-1.0
else
	B_m=B_m0
	!绝热膨胀时磁场减小，不考虑绝热膨胀则磁场为常数
endif

!由于本程序没有考虑空间扩散，我们做了假设，认为每一步一开始密度立即变成新时间的密度

!	if (j>1) then
!		uu=uu*((c_speed*(time-dt)*gamma_bulk+r0)/rr)**3
!	endif	

loop1: do i=1,n-1                        !注意这里的n,如果想算deltax的范围，要满足dx*n=deltax
!			print*,delta(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*********Note that we only solve points between 2-n. To 1 and n+1, N is given by the 
!boundary conditions!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (j==1) then
		x(i)=dx*i+loggamma_low	
	endif	
	
		Bi05=B(dx*(1.d0*i+0.5)+loggamma_low)    !B(i+0.5)
		Ci05=C(dx*(1.d0*i+0.5)+loggamma_low)    !C(i+0.5)
		Bi_05=B(dx*(1.d0*i-0.5)+loggamma_low)   !B(i-0.5)
		Ci_05=C(dx*(1.d0*i-0.5)+loggamma_low)   !C(i-0.5)

		
		!if (i==n+1) then
			!Bi05=B(dx*(1.d0*n)+loggamma_low)    !B(i+0.5)
			!Ci05=C(dx*(1.d0*n)+loggamma_low)    !C(i+0.5)
		!else  (i==1) then
			!Bi_05=B(loggamma_low)   !B(i-0.5)
			!Ci_05=C(loggamma_low)   !C(i-0.5)
		!else
			!Bi05=B(dx*(1.d0*(i-1)+0.5)+loggamma_low)    !B(i+0.5)
			!Ci05=C(dx*(1.d0*(i-1)+0.5)+loggamma_low)    !C(i+0.5)
			!Bi_05=B(dx*(1.d0*(i-1)-0.5)+loggamma_low)   !B(i-0.5)
			!Ci_05=C(dx*(1.d0*(i-1)-0.5)+loggamma_low)   !C(i-0.5)
		!endif

		w05=dx*Bi05/Ci05           !w(i+0.5)
		w_05=dx*Bi_05/Ci_05        !w(i-0.5)
	if 	(w05>100.) then
		delta05=1./w05
	else	
		delta05=1./w05-1./(dexp(w05)-1.d0)    !delta(i+0.5) 
	endif	
	
	if 	(w_05>100.) then
		delta_05=1./w_05
	else	
		delta_05=1./w_05-1./(dexp(w_05)-1.d0) !delta(i-0.5)
	endif

		C1(i)=(1.-delta05)*Bi05+1./dx*Ci05
	        C2(i)=-1./dx*(Ci05+Ci_05)-(1.d0-delta_05)*Bi_05+delta05*Bi05
		C3(i)=1.d0/dx*Ci_05-delta_05*Bi_05
!		print*,1.d0/dx*Ci_05-delta_05*Bi_05,1.d0/dx*Ci_05,delta_05*Bi_05
		A_A(i)=A(dx*i+loggamma_low)   !!!C1,C2,C3是系数，与note里一致。这些系数用到了Chang & Cooper文中的系数A(A_A),B,C等函数
				 !!!还用到delta函数，这个函数也用到A,B,C等三个函数。	

			if	(i>1) then
				aa(i)=C3(i)/A_A(i)
				if (abs(aa(i))<1.d-18) then
					aa(i)=0.d0
				endif
			endif
			
			if	(i<n-1) then
				cc(i)=C1(i)/A_A(i)
			endif
			
			bb(i)=C2(i)/A_A(i)-ratio       !AA, BB, CC, r与Press1992程序里三对角矩阵一致。
		if (j==1)		then
			
!				r(i)=-ratio*u(dexp(dx*i+loggamma_low))-dx*Q(dx*i+loggamma_low)           !u是初值函数, r付初始值
				r(i)=-dx*Q(dx*i+loggamma_low)                                             
		elseif (j<=int(inj_t/dt) ) then           !注入1/dt时间后停止
				r(i)=-ratio*uu(i)-dx*Q(dx*i+loggamma_low)*kk   !第二次计算起用第一次算的u来算r.
!				print*,ratio*uu(i),dx*Q(dx*i+loggamma_low)
		else
				r(i)=-ratio*uu(i)	!不再注入
		endif
!if (i==10) then
!endif
end do loop1
print*,j


!open(10001,file='r_v2.txt')
!do i=1,n-1
	!write(10001,*)log10(dexp(dx*i+loggamma_low)),r(i)
!enddo

call tridag(aa,bb,cc,r,uu,n-1)           !初始输入r,解出来uu即所要的结果


continue
end do loop2


loop3: do i=1,n-1
		Inuu=Inu(1.d0*10.**(10.+13.*i/(n-1.)))
!			write(*,*)10.+13.*i/(n-1.),dlog10(Inuu)
!			if (10.+13.*i/(n-1.)>20.5) then
!				continue
!			endif
		if (Inuu>1.d-20) then	
			write(4,*)10.+13.*i/(n-1.)+dlog10(gamma_bulk),dlog10(Inuu)
		endif	
		if (uu(i)>1.d-40) then
			write(*,*)log10(dexp(dx*i+loggamma_low)),log10(uu(i))
			write(3,*)log10(dexp(dx*i+loggamma_low)),log10(uu(i))
		endif	
	   enddo  loop3
end


!end of the main program


!--------------------------------------------------------------------------------------



recursive function Q(x)
	use initial_parameter
        implicit none
        real*8 x,Q,gammae,q0,rr
	common /qq0/  q0
	common /rr0/ rr
	gammae=dexp(x)
	if (gammae>gamma_low_inj.and.gammae<gamma_up_inj) then
		q=q0*gammae**(-p) !/(4.*3.1415926*rr**3*gamma_bulk**2)
	!	print*,(4.*3.1415926*rr**3*gamma_bulk**2)

	else 
		q=0.d0
	endif
end


recursive function u(gammae)             !电子谱初,值gamma_e
use initial_parameter
implicit none
real*8 gammae,u,dt,q
common /dtt/ dt
external q
if (gammae>gamma_low_inj.and.gammae<gamma_up_inj) then
	u=q(dlog(gammae))*dt
else 
	u=0.d0
endif
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function  A(x),B(x),C(x)   !The C below equation (7) of note.tex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive function A(x)           !The A below equation (7) of note.tex

implicit none
real*8 x,A
A=dexp(x)

end

recursive function B(x)           !The B below equation (7) of note.tex
        implicit none
        real*8 x,Csyn,CIC,H,B,Htmp
	real*8 Csynn,tmp,adiab,ab,ex
	external adiab
        external Csyn,CIC,H
	common /htmp/ Htmp
	common /abb/ ab
	Htmp=H(x) 		!与C函数共用一个H	
!       B=Csyn(x)+CIC(x)-2.*Htmp*dexp(-x)
	Csynn=Csyn(x)
!	ex=dexp(x)
!	tmp=-2.*Htmp*dexp(-x)
!        tmp=Htmp*(1.-2.*ex**2)/ex/(ex**2-1.)  
	B=Csynn+CIC(x)+ab*adiab(x)+tmp

!print*,Csynn,tmp
        end



recursive function adiab(x)
	use initial_parameter
	use constant
	implicit none
	real*8 x,adiab,r0,time
	common /tim/ time,r0
	adiab=dexp(x)/(time+r0/(c_speed*gamma_bulk))
	end	



recursive function C(x)             !The C below equation (7) of note.tex
implicit none
real*8 x,C,Htmp
!external H
common /htmp/ Htmp
	C=Htmp*dexp(-x)             !Htmp与B共用H（x）
end






recursive function Csyn(x)
use constant
use initial_parameter  !由主程序给出值
implicit none
real*8 x,gamma_e,Csyn,CICtmp
common /ccic/ CICtmp
gamma_e=dexp(x)

Csyn=4./3.*sigma_T*c_speed*(gamma_e*gamma_e-1.)*B_m**2/8./pi/m_e/c_speed**2
CICtmp=Csyn
end


recursive function CIC(x)
	use comptonY
	use constant
	use gaus_int
        implicit none
	integer i,icic
        real*8 x,CIC,CICtmp,Inu,nutmp,CIC_int(nnn),gamma_e,sigma
	real*8 t_tot,Uph
common /ccic/ CICtmp
external Inu,sigma
gamma_e=dexp(x)

!CIC=comp_Y*CICtmp  ! if we use a simplicity Y
!dndnu=Inu/(h_planck*nu)*4.*pi/c_speed !if use the precise form 
icic=0
!if (icic==1) then
	!do i=1,nnn
		!nutmp=nuuu(i)
                !CIC_int(i)=Inu(nutmp)*4.*pi/c_speed*c_speed*sigma(gamma_e,nutmp)
!!		print*,nnn,i,CIC_int(i)
        !enddo


	!CIC=sum(CIC_int*ww2)*4./3.*gamma_e**2/m_e/c_speed**2
!!	print*,dlog10(sum(CIC_int*ww2)/c_speed/sigma_T)
!else
        !CIC=comp_Y*CICtmp 

	!endif
	call Uph_gamma_c(t_tot,Uph)
CIC=4./3.*sigma_T*c_speed*gamma_e*gamma_e*Uph/m_e/c_speed**2


end

recursive subroutine Uph_gamma_c(t_tot,Uph) !solve equation (1) and (2) for gc and Uph
	use constant
	use initial_parameter
	implicit none
	real*8 Q,ti,t_tot,gamma_c,Uph,tau
	real*8 Uph_f,tol,zbrent,UB
	external Uph_f,zbrent,Q
	tol=1.d-7
!Uph=zbrent(Uph_f,0.d0,1.d16,tol) 
!UB=B_m**2/8./pi
Uph=1.d0
!gamma_c=3.*pi*m_e*c_speed/(4.*sigma_T*UB*(t_tot-ti)*dexp(-tau)+4.*sigma_T*Uph*(t_tot-ti)+3.*m_e*c_speed/gamma_up_inj)
	end

!recursive function Uph_f(Uph)
	!use constant
	!use initial_parameter
!implicit none
!integer,parameter::nnn=65
!real*8 Uph_f,t_tot,hat_t
!real*8 tt1(nnn),ww1(nnn),tmp(nnn)
!integer i
!real*8 gcc,int_hat_t,Uph,Q0,UB,tau,tt
!external tau
	!common /qq0/  q0
	!common /times/ t_tot


!gcc(hat_t)=3.*pi*m_e*c_speed/(4.*sigma_T*UB*(t_tot-hat_t)*dexp(-tau)+4.*sigma_T*Uph*(t_tot-hat_t)+3.*m_e&
!&*c_speed/gamma_up_inj)       !t_tot由主程序给定.最后的gamma_max是因为从gamma_max冷却到gamma_c与从2gamma_c冷却到gamma_c
                           !!的时间几乎相同   
!UB=B_m**2/8./pi

!call gaus_point(0.d0,t_tot,tt1,ww1,nnn) !tt1为时间节点
	!do i=1,nnn
		!tt=tt1(i)
		!tmp(i)=Q0*(gamma_up_inj**(2.-p)-gcc(tt)**(2.-p))/(2.-p) !（1）式第二个积分,p由主程序给定
        !enddo
!int_hat_t=g*m_e*c_speed**2*sum(tmp*ww1)
!Uph_f=int_hat_t-Uph
!!Uph_f=(3.*pi*m_e*c_speed/gamma_c-3.*m_e*c_speed/gamma_e-4.*sigma_T*UB*(t_tot-hat_t)*dexp(-tau))/4./sigma_T&
!!&/(t_tot-hat_t)-g*int_hat_t
	!end


!recursive function tau(gc)                    !辐射强度
!use initial_parameter  !由主程序给出值Radius和B_m
!use gaus_int           !integral points used
!use constant
!implicit none
!real*8 tau,gc
!!Real*8 jnu,alphanu,nu,absor
!real*8 nu,gPN,Npnu,alphanu_tmp(nnn)
!real*8 alphanu_s,alphanu
!integer i
!!external jnu,alphanu
!external gPn,Npnu

	!nu=q_e*B_m/(2.*pi*m_e*c_speed)*gc**2
	!do i=1,nnn
                !alphanu_tmp(i)=gPN(nu,xxx(i))
	!enddo

        !alphanu_s=sum(alphanu_tmp*ww1)
	!alphanu=1./(8.*pi*m_e*nu*nu*nu*h_planck/(m_e*c_speed**2))*alphanu_s
	
	!if (alphanu*Radius<=1.d-13) then
		!tau=0.d0
		!return
	!else if (alphanu*Radius>=100.) then 
                !tau=100.
		!return
	!else
		!tau=alphanu*Radius
	!endif	

!end


recursive function sigma(gamma_e,nutmp)
	use constant
	implicit none
real*8 gamma_e,nutmp,sigma,x

if (x<=0.000001) then
	sigma=sigma_T
else	
	x=gamma_e*h_planck*nutmp/m_e/c_speed**2

	sigma=3./4.*sigma_T*((1.+x)/x**3*(2.*x*(1.+x)/(1.+2.*x)-dlog(1.d0+2.*x))+1./(2.*x)&
	&*dlog(1.d0+2.*x)-(1.+3.*x)/(1.+2.*x)**2)
endif
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function H(x),Inu(nu),jnu(nu),alpha(nu)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive function H(x)
!use Hgamma_e
!use nu_up_low
use constant
use gaus_int
implicit none
real*8 I_Pnu_tmp(nnn),I_Pnu_s,resul,gamma_e,x,H,Hga_e
real*8 Inu,Pnu,nutmp,ggam,csyn
integer i,j
external Inu,Pnu,Csyn
Hga_e=dexp(x)

!open(123,file='test1.txt')
!open(124,file='test2.txt')

!do j=1,100
	!ggam=10.**(log10(3.2)+(log10(1.d4)-log10(3.2))*j/100.)
	!Hga_e=ggam

!print*,nu_low,nu_up

	do i=1,nnn
		nutmp=nuuu(i)
                I_Pnu_tmp(i)=Inu(nutmp)*Pnu(nutmp,Hga_e)/(2.*m_e*nutmp**2)
        enddo

        I_Pnu_s=sum(I_pnu_tmp*ww2)/(m_e*c_speed**2)

!call gaus(I_Pnu,nu_low,nu_up,resul,5)
H=I_Pnu_s

!if (H>=0.d0) then
	!write(123,*)dlog10(ggam),dlog10(H/ggam)
	!write(124,*)dlog10(ggam),dlog10(Csyn(dlog(ggam)))

!endif


!enddo
!close(123)
!close(124)
end



recursive function pnu(nu,gamma_e)         !谱功率，需要给定电子能量
use constant
use initial_parameter
implicit none
real*8 gamma_e,UB,nu_B,nu,nu_L,pnu
real*8 y
real*8 ri,rk1,rip,rkp,rk2
!nu_B=q_e*B_m/(2.*pi*m_e*c_speed)
nu_L=1./(2.*pi)*q_e*B_m/m_e/c_speed
y=nu/(3.*gamma_e**2*nu_L)
UB=B_m*B_m/8./pi
if (y<=60.d0) then
 call bessik (y,4.0d0/3.d0,ri,rk1,rip,rkp)                 !K贝塞尔函数是个减函数，到y=60时已经减小到1d-27
ri=0.d0; rip=0.d0; rkp=0.d0
 call bessik (y,1.0d0/3.0d0,ri,rk2,rip,rkp)
ri=0.d0; rip=0.d0; rkp=0.d0 
else
rk1=0.d0
rk2=0.d0
endif
!Pnu=3.*3.**0.5/pi*sigma_T*c*UB*y**2*( Kb(4./3.,y)*Kb(1./3.,y)-3./5.*y*(Kb(4./3.,y)**2-Kb(1./3.,y)**2) )
!Pnu=3.*3.**0.5/pi*sigma_T*c*UB*y**2*( rk1*rk2-3./5.*y*(rk1**2-rk2**2) ) !Chen 2011
Pnu=4.*pi*3.**0.5*q_e**2*nu_L/c_speed*y**2*( rk1*rk2-3./5.*y*(rk1**2-rk2**2) )
end


recursive function Inu(nu)                    !辐射强度
use initial_parameter  !由主程序给出值Radius和B_m
use gaus_int           !integral points used
use constant
implicit none
real*8 Inu
!Real*8 jnu,alphanu,nu,absor
real*8 nu,gPN,Npnu,alphanu_tmp(nnn),jnu_tmp(nnn)
real*8 alphanu_s,alphanu,jnu_s,jnu
integer i,kabs
!external jnu,alphanu
external gPn,Npnu
common /kabss/ kabs
!open(123,file='test.txt')
!do j=1,100
	!nu=10.**(9.+(24.-9.)*j/100.)

	do i=1,nnn
                alphanu_tmp(i)=gPN(nu,xxx(i))
		    jnu_tmp(i)=Npnu(nu,xxx(i))
        enddo

        alphanu_s=sum(alphanu_tmp*ww1)
	alphanu=1./(8.*pi*m_e*nu*nu*nu*h_planck/(m_e*c_speed**2))*alphanu_s
	jnu_s=sum(jnu_tmp*ww1)
	jnu=1./4./pi*jnu_s
	
!print*,1.d0-dexp(-1.d-13)
!absor=alphanu(nu)
!if (alphanu_s<=0.d0) then
!	print*,"something wrong"
!	Inu=jnu*Radius
!	Inu=0.d0
!else 
!	print*,alphanu*radiuhi
	if (alphanu*Radius<=1.d-13) then
		Inu=jnu*Radius
	elseif (alphanu*Radius>100.) then         !指数不能太大
!			print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!			print*,alphanu*Radius,dexp(-alphanu*Radius)
!			pause
!		endif
			Inu=jnu/alphanu           
	else
		Inu=jnu/alphanu*(1.-dexp(-alphanu*Radius))

	endif

	if (kabs==0) Inu=jnu*Radius    !if kabs==0对应不考虑自吸收
		!	print*,jnu*Radius,Inu
!endif

!if (Inu>=0.d0) then
!	write(123,*)dlog10(nu),dlog10(Inu)
!endif

!enddo


end




recursive function Npnu(nu,gamma_e)    !发射系数的被积函数N*Pnu
!use jnu_nu
implicit none
real*8 gamma_e,Nga,pnu,u,Npnu,nu,utmp
integer j
common /j/ j
external Nga,pnu,u

if (j==1) then
	utmp=u(gamma_e)
	if (utmp>=0.d0) then
		Npnu=utmp*pnu(nu,gamma_e)      !初始时刻通过初值得到 
	else
		Npnu=0.d0
	endif
else	
	Npnu=Nga(gamma_e)*pnu(nu,gamma_e)    !Nga通过插值得到任意一点的电子谱,nu通过jnu_nu从发射系数里得到
endif
end



recursive function gPN(nu,gamma_e)  !吸收系数中的被积函数
!use jnu_nu
use constant
use initial_parameter
implicit none
real*8 gamma_e,Pnu,h_in,err,E_e,Ntmp,Nga,u,gPN,nu
real*8 deviative_func,dfridr
integer j
external Pnu,Nga,u,dfridr,deviative_func
common /j/ j
!Note if we use dfridr, then it will be more precise, but more consuming time.
!If we use the current form, note that in Inu, the nu's power is 3. If we use 
!dfridr, then it is 2.

if (j==1) then

	!Ntmp=u(gamma_e)/gamma_e**2    !需要用到某一时刻电子分布函数,初始时刻是个函数，后来是计算出的结果
	if (gamma_e-h_planck*nu/m_e/c_speed**2>=gamma_low) then
		!光子能量必须小于电子能量
		gPN=gamma_e*gamma_e*pnu(nu,gamma_e)*(u(gamma_e-h_planck*nu/m_e/c_speed**2)/&
		&(gamma_e-h_planck*nu/m_e/c_speed**2)**2-u(gamma_e)/gamma_e**2)

       !         if (gPN<0.d0) then
			!print*,gamma_e-h_planck*nu/m_e/c_speed**2,u(gamma_e-h_planck*nu/m_e/c_speed**2)/&
		!&(gamma_e-h_planck*nu/m_e/c_speed**2)**2,u(gamma_e)/gamma_e**2
			!continue
		!endif
	else
		gPN=0.d0
	endif	

else

	!Ntmp=Nga(gamma_e)/gamma_e**2   
	if (gamma_e-h_planck*nu/m_e/c_speed**2>=gamma_low) then
		gPn=gamma_e*gamma_e*pnu(nu,gamma_e)*(Nga(gamma_e-h_planck*nu/m_e/c_speed**2)/&
		&(gamma_e-h_planck*nu/m_e/c_speed**2)**2-Nga(gamma_e)/gamma_e**2)
 !       	if ((Nga(gamma_e-h_planck*nu/m_e/c_speed**2)/(gamma_e-h_planck*nu/m_e/c_speed**2)**2-Nga(gamma_e)/gamma_e**2)<0.d0) then
!			print*,"pause"
!			print*,gamma_e
!			print*,(Nga(gamma_e-h_planck*nu/m_e/c_speed**2)/(gamma_e-h_planck*nu/m_e/c_speed**2)**2-Nga(gamma_e)/gamma_e**2),&
!			&dfridr(deviative_func,gamma_e,0.01d0,err)*(-h_planck*nu)/m_e/c_speed**2 
!			print*,"pause stop"
!			pause
!		endif
!		if (gPN<0.d0) then
!			print*,gamma_e-h_planck*nu/m_e/c_speed**2,Nga(gamma_e-h_planck*nu/m_e/c_speed**2)/&
!			&(gamma_e-h_planck*nu/m_e/c_speed**2)**2,Nga(gamma_e)/gamma_e**2
!			continue
!		endif
	
	
	else
		gPN=0.d0
	endif

 !解出来的是固定点，需要用到插值得到任一点的值,in是电子分布函数的点数，要和主程序保持一致
endif	!gamma_ep是电子能量矩阵, Ngp是电子分布矩阵，要与主程序保持一致		

if (gPN<0.d0) gPN=0.d0

end            																		   



real*8 function deviative_func(gamma_e)             !微分的函数
implicit none
real*8 u,Nga,gamma_e
integer j
external u,Nga
common /j/ j                                        !初始时刻用初值

if (j==1) then
	deviative_func=u(gamma_e)/gamma_e**2                !需要用到某一时刻电子分布函数,初始时刻是个函数，后来是计算出的结果
else
	deviative_func=Nga(gamma_e)/gamma_e**2              !解出来的是固定点，需要用到插值得到任一点的值,in是电子分布函数的点数，要和主程序保持一致
endif												!gamma_ep是电子能量矩阵, Ngp是电子分布矩阵，要与主程序保持一致		


end




recursive function Nga(gamma_e) 			!interpolation
use Ng
implicit none 
integer nn
real*8 gamma_e,y,Nga,loggam
nn=1
!xx=dexp(x)
loggam=dlog(gamma_e)
y=0.d0
call Interpol_Akima(n-1,nn,loggam,y,x,uu)
Nga=y
end





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



real*8 FUNCTION dfridr(func,x,h,err)                !用理查德森外推法求微分
implicit none
      INTEGER NTAB
      REAL*8 err,h,x,func,CON,CON2,BIG,SAFE
      PARAMETER (CON=1.4,CON2=CON*CON,BIG=1.d30,NTAB=10,SAFE=2.)
      EXTERNAL func
!CU    USES func
      INTEGER i,j
      REAL*8 errt,fac,hh,a(NTAB,NTAB)
      if(h.eq.0.) pause 'h must be nonzero in dfridr'
      hh=h
      a(1,1)=(func(x+hh)-func(x-hh))/(2.0*hh)
      err=BIG
      do 12 i=2,NTAB
        hh=hh/CON
        a(1,i)=(func(x+hh)-func(x-hh))/(2.0*hh)
        fac=CON2
        do 11 j=2,i
          a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.)

          fac=CON2*fac
          errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
          if (errt.le.err) then
            err=errt
            dfridr=a(j,i)
          endif
11      continue
        if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err)return
12    continue
      return
      END







recursive subroutine gaus_point(x0,x2,y,w,n)          !高斯节点
implicit none
real*8 x_11,x_12,x_13,x_14,x_15,x_16,x_17,x_18,x_19,x_110,x_111,x_112,x_113,x1,x2,x0,step
integer nn,j,n

real*8 w(n),y(n)

external func
  
   if (x0==0.d0) then
   x1=1.d-9
   else
   x1=x0
   endif

   nn=n/13             

    step=log10(x2/x1)/13.
   	x_11=x1
    x_12=10.**(log10(x1)+1.*step)
	x_13=10.**(log10(x1)+2.*step)
    x_14=10.**(log10(x1)+3.*step)
	x_15=10.**(log10(x1)+4.*step)
	x_16=10.**(log10(x1)+5.*step)
	x_17=10.**(log10(x1)+6.*step)
	x_18=10.**(log10(x1)+7.*step)
	x_19=10.**(log10(x1)+8.*step)
	x_110=10.**(log10(x1)+9.*step)
	x_111=10.**(log10(x1)+10.*step)
	x_112=10.**(log10(x1)+11.*step)
	x_113=10.**(log10(x1)+12.*step)
 	call gauleg(x_11,x_12,y(1:nn),w(1:nn),nn)
	call gauleg(x_12,x_13,y(nn+1:2*nn),w(nn+1:2*nn),nn)
	call gauleg(x_13,x_14,y(2*nn+1:3*nn),w(2*nn+1:3*nn),nn)
	call gauleg(x_14,x_15,y(3*nn+1:4*nn),w(3*nn+1:4*nn),nn)
	call gauleg(x_15,x_16,y(4*nn+1:5*nn),w(4*nn+1:5*nn),nn)
	call gauleg(x_16,x_17,y(5*nn+1:6*nn),w(5*nn+1:6*nn),nn)
	call gauleg(x_17,x_18,y(6*nn+1:7*nn),w(6*nn+1:7*nn),nn)
	call gauleg(x_18,x_19,y(7*nn+1:8*nn),w(7*nn+1:8*nn),nn)
	call gauleg(x_19,x_110,y(8*nn+1:9*nn),w(8*nn+1:9*nn),nn)
        call gauleg(x_110,x_111,y(9*nn+1:10*nn),w(9*nn+1:10*nn),nn)
        call gauleg(x_111,x_112,y(10*nn+1:11*nn),w(10*nn+1:11*nn),nn)
	call gauleg(x_112,x_113,y(11*nn+1:12*nn),w(11*nn+1:12*nn),nn)
	call gauleg(x_113,x2,y(12*nn+1:13*nn),w(12*nn+1:13*nn),nn)
!resul=0.
!do j=1,nn
!     resul=w1(j)*func(y1(j))+w2(j)*func(y2(j))+w3(j)*func(y3(j))+&
!&	 w4(j)*func(y4(j))+w5(j)*func(y5(j))+w6(j)*func(y6(j))+w7(j)*func(y7(j))+&
!&    w8(j)*func(y8(j))+w9(j)*func(y9(j))+w10(j)*func(y10(j))+w11(j)*func(y11(j))+&
!&    w12(j)*func(y12(j))+w13(j)*func(y13(j))+resul
!   enddo

  !  deallocate(w1,y1,w2,y2,w3,y3,w4,y4,w5,y5,w6,y6,w7,y7)
!  deallocate(w8,y8,w9,y9,w10,y10))
  end






recursive subroutine gaus(func,x0,x2,resul,nn)
implicit none
real*8 x_11,x_12,x_13,x_14,x_15,x_16,x_17,x_18,x_19,x_110,x_111,x_112,x_113,x1,x2,resul1,func,x0
integer nn,j
real*8 w1(nn),y1(nn),w2(nn),y2(nn),w3(nn),y3(nn),w4(nn),y4(nn),w5(nn),y5(nn),w6(nn),y6(nn),w7(nn),y7(nn),step,resul
real*8 w8(nn),y8(nn),w9(nn),y9(nn),w10(nn),y10(nn),y11(nn),y12(nn),y13(nn),w11(nn),w12(nn),w13(nn)
external func
  
   if (x0==0.d0) then
   x1=1.d-9
   else
   x1=x0
   endif

    step=log10(x2/x1)/13.
   	x_11=x1
    x_12=10.**(log10(x1)+1.*step)
	x_13=10.**(log10(x1)+2.*step)
    x_14=10.**(log10(x1)+3.*step)
	x_15=10.**(log10(x1)+4.*step)
	x_16=10.**(log10(x1)+5.*step)
	x_17=10.**(log10(x1)+6.*step)
	x_18=10.**(log10(x1)+7.*step)
	x_19=10.**(log10(x1)+8.*step)
	x_110=10.**(log10(x1)+9.*step)
	x_111=10.**(log10(x1)+10.*step)
	x_112=10.**(log10(x1)+11.*step)
	x_113=10.**(log10(x1)+12.*step)
 	call gauleg(x_11,x_12,y1,w1,nn)
	call gauleg(x_12,x_13,y2,w2,nn)
	call gauleg(x_13,x_14,y3,w3,nn)
	call gauleg(x_14,x_15,y4,w4,nn)
	call gauleg(x_15,x_16,y5,w5,nn)
	call gauleg(x_16,x_17,y6,w6,nn)
	call gauleg(x_17,x_18,y7,w7,nn)
	call gauleg(x_18,x_19,y8,w8,nn)
	call gauleg(x_19,x_110,y9,w9,nn)
    call gauleg(x_110,x_111,y10,w10,nn)
    call gauleg(x_111,x_112,y11,w11,nn)
	call gauleg(x_112,x_113,y12,w12,nn)
	call gauleg(x_113,x2,y13,w13,nn)
resul=0.
do j=1,nn
     resul=w1(j)*func(y1(j))+w2(j)*func(y2(j))+w3(j)*func(y3(j))+&
&	 w4(j)*func(y4(j))+w5(j)*func(y5(j))+w6(j)*func(y6(j))+w7(j)*func(y7(j))+&
&    w8(j)*func(y8(j))+w9(j)*func(y9(j))+w10(j)*func(y10(j))+w11(j)*func(y11(j))+&
&    w12(j)*func(y12(j))+w13(j)*func(y13(j))+resul
   enddo

  !  deallocate(w1,y1,w2,y2,w3,y3,w4,y4,w5,y5,w6,y6,w7,y7)
!  deallocate(w8,y8,w9,y9,w10,y10))
  end
recursive SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 112 i=1,m
        z=dcos(3.141592654d0*(i-.25d0)/(n+.5d0))
1000       continue
          p1=1.d0
          p2=0.d0
          do 111 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
111        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z

          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1000
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
112    continue
      return
      END




	  recursive  SUBROUTINE tridag(a,b,c,r,u,n)
	  implicit none
      INTEGER n,NMAX
      REAL*8 a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=90000)
      INTEGER j
      REAL*8 bet,gam(NMAX)
      if(b(1).eq.0.)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 113 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
113    continue
      do 123 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
123    continue
      return
      END







	  
!********************************************************
!*          Akima spline fitting subroutine             *
!* ---------------------------------------------------- *
!* The input table is X(i), Y(i), where Y(i) is the     *
!* dependant variable. The interpolation point is xx,   *
!* which is assumed to be in the interval of the table  *
!* with at least one table value to the left, and three *
!* to the right. The interpolated returned value is yy. *
!* n is returned as an error check (n=0 implies error). *
!* It is also assumed that the X(i) are in ascending    *
!* order.                                               *
!********************************************************
recursive Subroutine Interpol_Akima(iv,n,xx,yy,X,Y)  
implicit none
  !Labels: 100,200,300
  integer size
  parameter(SIZE=9999)
  integer i,iv,n
  real*8  xx,yy
  real*8  X (0:SIZE), Y (0:SIZE)
  real*8  XM(0:SIZE+3)
  real*8  Z (0:SIZE)
  real*8 a,b

  n=1
  !special case xx=0
  if (xx.eq.0.0) then
    yy=0.d0; return
  end if
  !Check to see if interpolation point is correct
  if (xx<X(1).or.xx>=X(iv-3)) then
    n=0 ; return
  end if
  X(0)=2.d0*X(1)-X(2)
  !Calculate Akima coefficients, a and b
  do i = 1, iv-1
    !Shift i to i+2
    XM(i+2)=(Y(i+1)-Y(i))/(X(i+1)-X(i))
  end do
  XM(iv+2)=2.d0*XM(iv+1)-XM(iv)
  XM(iv+3)=2.d0*XM(iv+2)-XM(iv+1)
  XM(2)=2.d0*XM(3)-XM(4)
  XM(1)=2.d0*XM(2)-XM(3)
  do i = 1, iv
    a=dabs(XM(i+3)-XM(i+2))
    b=dabs(XM(i+1)-XM(i))
    if (a+b.ne.0.d0) goto 100
    Z(i)=(XM(i+2)+XM(i+1))/2.d0
    goto 200
100 Z(i)=(a*XM(i+1)+b*XM(i+2))/(a+b)
200 end do
  !Find relevant table interval
  i=0
300 i=i+1
  if (xx>X(i)) goto 300
  i=i-1
  !Begin interpolation
  b=X(i+1)-X(i)
  a=xx-X(i)
  yy=Y(i)+Z(i)*a+(3.d0*XM(i+2)-2.d0*Z(i)-Z(i+1))*a*a/b
  yy=yy+(Z(i)+Z(i+1)-2.d0*XM(i+2))*a*a*a/(b*b)
  return
end


recursive SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
!*****************************************************
!*     Polynomial Interpolation or Extrapolation     *
!*            of a Discreet Function                 *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*   XA:    Table of abcissas  (N)                   *
!*   YA:    Table of ordinates (N)                   *
!*    N:    Number of points                         *
!*    X:    Interpolating abscissa value             *
!* OUTPUTS:                                          *
!*    Y:    Returned estimation of function for X    *
!*   DY:    Estimated error for Y                    *
!*****************************************************
implicit none
integer nmax,n,ns,i,m
PARAMETER(NMAX=200)
REAL*8 XA(N),YA(N), X,Y,DY
REAL*8 C(NMAX),D(NMAX)
REAL*8 DEN,DIF,DIFT,HO,HP,W
NS=1
DIF=DABS(X-XA(1))
DO I=1,N
  DIFT=DABS(X-XA(I))
  IF(DIFT.LT.DIF) THEN
    NS=I                 !index of closest table entry
	DIF=DIFT
  ENDIF
  C(I)=YA(I)             !Initialize the C's and D's
  D(I)=YA(I)
END DO
Y=YA(NS)                 !Initial approximation of Y
NS=NS-1
DO M=1,N-1
  DO I=1,N-M
    HO=XA(I)-X
	HP=XA(I+M)-X
    W=C(I+1)-D(I) 
    DEN=HO-HP
	IF(DEN.EQ.0.) Pause 'Error: two identical abcissas)'
	DEN=W/DEN
	D(I)=HP*DEN          !Update the C's and D's
	C(I)=HO*DEN
  END DO
  IF(2*NS.LT.N-M) THEN   !After each column in the tableau XA
    DY=C(NS+1)           !is completed, we decide which correction,
  ELSE                   !C or D, we want to add to our accumulating 
    DY=D(NS)             !value of Y, i.e. which path to take through 
	NS=NS-1              !the tableau, forking up or down. We do this
  ENDIF                  !in such a way as to take the most "straight 
  Y=Y+DY	             !line" route through the tableau to its apex,
END DO                   !updating NS accordingly to keep track of 
                         !where we are. This route keeps the partial
RETURN                   !approximations centered (insofar as possible)
END                      !on the target X.The last DY added is thus the
                         !error indication.

!end of file tpolint.f90








recursive FUNCTION chebev(a,b,c,m,x)
      implicit none
      INTEGER m
      real*8 chebev,a,b,x,c(m)
      INTEGER j
      real*8 d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.d0) pause 'x not in range in chebev'
      d=0.
      dd=0.
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do 18 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
18    continue
      chebev=y*d-dd+0.5d0*c(1)
      return
      END




recursive SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      implicit none
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=5,NUSE2=5)
      !USE chebev
      real*8 xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371168d0,6.5165112670737d-3,3.087090173086d-4,-3.4706269649d-6,6.9437664d-9,3.67795d-11,-1.356d-13/
      DATA c2/1.843740587300905d0,-7.68528408447867d-2,1.2719271366546d-3,-4.9717367042d-6,-3.31261198d-8,2.423096d-10,-1.702d-13,-1.49d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END


recursive SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
      implicit none
     !double precision,external::beschb
      INTEGER MAXIT
      real*8 ri,rip,rk,rkp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.e-6,FPMIN=1.e-10,MAXIT=10000,XMIN=2.d0,PI=3.141592653589793d0)
      !USE beschb
      INTEGER i,l,nl
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessik'
      nl=int(xnu+.5d0)
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=1.d0/(b+d)
        c=b+1.d0/c
        del=c*d
        h=del*h
        if(abs(del-1.d0).lt.EPS)goto 1
!		print*,del,x
11    continue
      pause 'x too large in bessik; try asymptotic expansion'
1     continue
      ril=FPMIN
      ripl=h*ril
      ril1=ril
      rip1=ripl
      fact=xnu*xi
      do 12 l=nl,1,-1
        ritemp=fact*ril+ripl
        fact=fact-xi
        ripl=fact*ritemp+ril
        ril=ritemp
12    continue
      f=ripl/ril
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=fact*(gam1*cosh(e)+gam2*fact2*d)
        sum=ff
        e=exp(e)
        p=0.5d0*e/gampl
        q=0.5d0/(e*gammi)
        c=1.d0
        d=x2*x2
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*ff
          sum=sum+del
          del1=c*(p-i*ff)
          sum1=sum1+del1
          if(abs(del).lt.abs(sum)*EPS)goto 2
13      continue
        pause 'bessk series failed to converge'
2       continue
        rkmu=sum
        rk1=sum1*xi2
      else
        b=2.d0*(1.d0+x)
        d=1.d0/b
        delh=d
        h=delh
        q1=0.d0
        q2=1.d0
        a1=.25d0-xmu2
        c=a1
        q=c
        a=-a1
        s=1.d0+q*delh
        do 14 i=2,MAXIT
          a=a-2*(i-1)
          c=-a*c/i
          qnew=(q1-b*q2)/a
          q1=q2
          q2=qnew
          q=q+c*qnew
          b=b+2.d0
          d=1.d0/(b+a*d)
          delh=(b*d-1.d0)*delh
          h=h+delh
          dels=q*delh
          s=s+dels
          if(abs(dels/s).lt.EPS)goto 3
14      continue
        pause 'bessik: failure to converge in cf2'
3       continue
        h=a1*h
        rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s
        rk1=rkmu*(xmu+x+.5d0-h)*xi
      endif
      rkmup=xmu*xi*rkmu-rk1
      rimu=xi/(f*rkmu-rkmup)
      ri=(rimu*ril1)/ril
      rip=(rimu*rip1)/ril
      do 15 i=1,nl
        rktemp=(xmu+i)*xi2*rk1+rkmu
        rkmu=rk1
        rk1=rktemp
15    continue
      rk=rkmu
      rkp=xnu*xi*rkmu-rk1
      return 
      END






     recursive FUNCTION zbrent(func,x1,x2,tol) !tol is the accuracy.
      INTEGER ITMAX
      REAL*8 zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL*8 a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause &
     'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa

          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5d0*tol
        xm=.5d0*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else

            q=fa/fc
            r=fb/fc
            p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.d0*p .lt. min(3.d0*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else

          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END

