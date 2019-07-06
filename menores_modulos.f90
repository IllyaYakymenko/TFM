! #################################################################################
! #
! #   Obtención de datos sobre los menores módulos
! #   
! #   Input: Channel 10 (fichero _raman.dat)
! #   Output: Channel 20: Valores de los menores módulos 
! #
! #################################################################################

program menores_modulos

implicit none

! Input

integer :: numfiles, nrow, ncol, numaa
character*200, allocatable :: filenames(:)
real :: perc

! Variables de análisis

real, allocatable :: diedros(:), diedrosant(:), desps(:)
double precision :: dPhi, dPsi, maxmodulo, modulodesp, cmaxmodulo
integer :: t

! Variables de iteración

integer :: fileind, rowiter, di, aa



read(5,*) numfiles
allocate(filenames(numfiles))
read(5,*) filenames
read(5,*) nrow
read(5,*) ncol
read(5,*) perc

numaa = ncol/2

allocate(diedros(ncol))
allocate(diedrosant(ncol))
allocate(desps(ncol))

maxmodulo = 0.d0

! Detectar máximo módulo (primer bucle)

print*, "Detectando módulo máximo"

do fileind = 1, numfiles
	open(10, file = filenames(fileind))
	read(10,*) diedrosant

	do rowiter=2, nrow
		read(10,*) diedros

		do di=1, ncol
			desps(di) = (diedros(di) - diedrosant(di))-1*int((diedros(di) - diedrosant(di))/180)*360
		end do

		do aa=1, numaa
			dPhi = desps(2*aa - 1)
			dPsi = desps(2*aa)

			modulodesp = sqrt(dPhi*dPhi + dPsi*dPsi)
			maxmodulo = max(modulodesp, maxmodulo)
		end do

		diedrosant = diedros

	end do
	print '("Detectando módulo mínimo: [", i3, "/", i3, "]")', fileind, numfiles

	close(10)
end do

! Establecimiento del límite de registro
! Se registran módulos cuyo valor sea <= a cmaxmodulo
cmaxmodulo = maxmodulo * (perc/100)

! Registro de valores cercanos a los mínimos (segundo bucle)

print*, "Registrando valores"

do fileind = 1, numfiles
	open(10, file = filenames(fileind))
	read(10,*) diedrosant

	t = 0

	do rowiter=2, nrow
		read(10,*) diedros

		t = t+1

		do di=1, ncol
			desps(di) = (diedros(di) - diedrosant(di))-1*int((diedros(di) - diedrosant(di))/180)*360
		end do

		do aa=1, numaa
			dPhi = desps(2*aa - 1)
			dPsi = desps(2*aa)

			modulodesp = sqrt(dPhi*dPhi + dPsi*dPsi)

			! Comprobación del tamaño del módulo para su escritura en el archivo de salida
			if(modulodesp.le.cmaxmodulo) then
				write(20,*) modulodesp, aa, diedrosant(2*aa -1), diedrosant(2*aa), diedros(2*aa-1), diedros(2*aa), dPhi, dPsi, t
			end if
			

		end do

		diedrosant = diedros

	end do
	print '("Registrando valores: [", i3, "/", i3, "]")', fileind, numfiles

	close(10)
end do

end program menores_modulos
