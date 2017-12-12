Program main

	Implicit none

	double precision :: tau, dt, T, R, tmax
	integer, parameter :: Np=500000
	double precision, parameter :: PI= 4*atan(1.)
	double precision, dimension(Np) :: Ep
	integer :: i,j
	!double precision, dimension(1000000) :: test
	double precision :: moy, start, finish, sigma
	double precision, dimension(:), allocatable :: Tint

	! Déclaration des paramètres
	Call CPU_TIME(start)
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

	open(unit=11, file="Tint.dat")
	Write(11,*) 0 , Tint(1)

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


		Write(11,*) i*dt, Tint(i+1)

	End do
	close(11)
	call histogramme("Histogramme",Ep)
	Call CPU_TIME(finish)
	Print*, "Temps de calcul : ", finish-start

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
		double precision ::  x, y, s
		double precision, intent(out) :: Z0

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
		double precision :: U1, U2
		double precision, intent(in):: PI
		double precision, intent(out) :: Z0

		call random_number(U1)
		call random_number(U2)
		Z0 = sqrt(-2*log(U1))*cos(2*PI*U2)


	end subroutine

	subroutine histogramme(name_file,Z)
		implicit none
		character(len=11), intent(in) :: name_file
		double precision, dimension(:), intent(in) :: Z
		integer :: i,j
		integer, dimension(100) :: C
		double precision  :: Emax, Emin, palier

		C=0
		Emin = Z(1)
		Emax = Z(1)
		Do i=1, size(Z)
			if (Z(i)>Emax) then
				Emax = Z(i)
			end if
			if (Z(i)<Emin) then
				Emin = Z(i)
			end if
		end do
		palier = (Emax - Emin)/100
		do i=1, size(Z)
			do j = 1, 100
				if ((Z(i) > Emin+(j-1)*palier) .and. (Z(i)<Emin+(j)*palier)) then
					C(j) = C(j) + 1
				end if
			end do
		end do

		open(unit=10, file=name_file)
		Do i=1,size(C)

			Write(10,*) (i-1)*palier+Emin, C(i)
			Write(10,*) (i-1)*palier+Emin+palier*0.99, C(i)
		end do
		close(10)


	end subroutine


End program main
