Program main

Implicit none

real :: tau, dt, T, R, time, tmax
integer, parameter :: Np=1000
real, parameter :: PI= 4*atan(1.)
real, dimension(Np) :: Ep
integer :: i,j
character(len=30) :: name_file
!real, dimension(1000000) :: test 
real :: var, moy, start, finish, sigma
real, dimension(:), allocatable :: Tint

! Déclaration des paramètres

tau = 1
dt = 0.01
T = 300
R = 280
tmax = 5

call random_number(Ep)
! Renvoie un nombre aléatoire entre 0 et 1

Ep = Ep*100000. + 100000.

moy = 0

Allocate(Tint(floor(tmax/dt)+1))

Do j=1,size(Ep)
	moy = moy + Ep(j)
End do
moy = moy/size(Ep)

Tint(1) = moy/R

Print*, "Tint à l'état initial = ", Tint(1), " moyenne = ", moy

Do i = 1,floor(tmax/dt)
	Do j = 1,size(Ep)
		call gaussienne2(sigma,PI)
		Ep(j) = (1/(1+2*dt/tau))*(Ep(j)+(R*T*dt/tau)*(1+sigma**2)+2*sqrt((dt/tau)*R*T*Ep(j))*sigma)
	End do
	moy = 0
	Do j=1,size(Ep)
		moy = moy + Ep(j)
	End do
	moy = moy/size(Ep)

	Tint(i+1) = moy/R
	Print*, "Tint(",i*dt,") = ", Tint(i+1)
End do

!-----------------------------------------------------------------------------------------------------------------------------

! Gaussienne 1

!Call CPU_TIME(start)

!do i=1,size(test)
!	call gaussienne1(test(i))
!end do

!Call CPU_TIME(finish)

!Print*, "Temps de la premiere methode : ", finish-start

!name_file = "gaussienne1.dat"
!call histogramme(name_file,test)

! Gaussienne 2

!Call CPU_TIME(start)

!do i=1,size(test)
!	call gaussienne2(test(i),PI)
!end do

!Call CPU_TIME(finish)

!Print*, "Temps de la deuxieme methode : ", finish-start

!name_file = "gaussienne2.dat"
!call histogramme(name_file,test)

!moy = 0
!var = 0

!Do j=1,size(test)
!	moy = moy + test(j)
!End do
!moy = moy/size(test)

!Do j=1,size(test)
!	var = var + (test(j)-moy)**2
!End do
!var = var/size(test)

!Print*, "Moyenne = ", moy
!Print*, "Variance = ", var



contains

subroutine gaussienne1(Z0)
	implicit none
	real ::  x, y, s
	real, intent(out) :: Z0

	call random_number(x)
	call random_number(y)
	x = x*2-1
	y = y*2-1
	s = x**2 + y**2
	do while ((s == 0) .or. (s >= 1)) 
		call random_number(x)
		call random_number(y)
		x = x*2-1
		y = y*2-1
		s = x**2 + y**2
	end do
	Z0 = x*sqrt((-2*log(s))/s)

end subroutine

subroutine gaussienne2(Z0, PI)
	implicit none
	real :: U1, U2
	real, intent(in):: PI
	real, intent(out) :: Z0

	call random_number(U1)
	call random_number(U2)
	Z0 = sqrt(-2*log(U1))*cos(2*PI*U2)


end subroutine

subroutine histogramme(name_file,Z)
	implicit none
	character(len=30), intent(in) :: name_file
	real, dimension(:), intent(in) :: Z
	integer :: i,j
	integer, dimension(1000) :: C
	real  :: p

	C=0
	
	Do i=1, size(Z)
		p = Z(i)*100
		j = floor(p)+501
		if ((j > 0) .and. (j<1000)) then 		
			C(j) = C(j)+1
		end if
	end do	

	open(unit=10, file=name_file)
	Do i=1,size(C)
		
		Write(10,*) i, C(i)
	end do
	close(10)
		

end subroutine


End program main
