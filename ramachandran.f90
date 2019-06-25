!	###############################################################
!	#
!	#   Obtención del diagrama de Ramachandran para cada
!	#   aminoácido.
!	#   Input: Channel 10 (ficheros _raman.dat)
!	#	Output: Channel 20: Ramachandran de todos los aminoácidos
!	#					21-> 20 + número de aminoácidos
!	#						Ramachandran de cada aminoácido
!	#
!	###############################################################


program ramachandran

implicit none

! Datos de entrada (archivo.in)

integer :: numfiles, resol, numrow, numcol
character*200, allocatable :: filenames(:)
real :: PhiMin, PhiMax, PsiMin, PsiMax

! Definición de la matriz del diagrama

double precision, allocatable :: diagrama(:,:,:)
real :: rangoPhi, rangoPsi
integer :: diagPhi, diagPsi, numaa
real, allocatable :: diedros(:)
integer, allocatable :: outputs(:)
double precision, allocatable :: totalpoints(:)

! Variables de control

integer ch, row, col, fileind, rowiter, diagrow, diagcol, aa, i, diagaa, waa

read(5,*) numfiles

allocate(filenames(numfiles))

read(5,*) filenames
read(5,*) numrow
read(5,*) numcol
read(5,*) PhiMin, PhiMax, PsiMin, PsiMax
read(5,*) resol

numaa = numcol/2

allocate(diagrama(numaa+1,resol,resol))
allocate(diedros(numcol))
allocate(outputs(numaa))
allocate(totalpoints(numaa+1))

! Creación de canales de output
! En outputs se guardan por orden las salidas de cada amonoácido
! Hay tantos outputs como aminoácidos en el polipéptido
! Las frecuencias de todos los aminoácidos se guardan en la posición
! numaa+1 de la primera dimensión de diagrama

do ch=1, numaa
	outputs(ch) = 20+ch	
end do

! Matriz del diagrama vacía
do diagaa=1, numaa+1
	do row=1, resol
		do col=1, resol
			diagrama(diagaa,row, col) = 0.d0
		end do
	end do
end do

rangoPhi = (PhiMax - PhiMin)/resol
rangoPsi = (PsiMax - PsiMin)/resol

do aa=1, numaa+1
	totalpoints(aa) = 0.d0
end do

print*, 'Calculando el diagrama de Ramachandran'

do fileind=1, numfiles
	open(10, file = filenames(fileind))
	do rowiter=1, numrow
		read(10,*)diedros
		do aa=1, numaa
			diagPhi = int((diedros(2*aa-1)+PhiMax)/rangoPhi) +1
			diagPsi = int((diedros(2*aa)+PsiMax)/rangoPsi)+1
			if(diagPhi.ge.resol) diagPhi=resol
			if(diagPsi.ge.resol) diagPsi=resol

			diagrama(aa, diagPsi, diagPhi) = diagrama(aa, diagPsi, diagPhi) +1
			diagrama(numaa+1, diagPsi, diagPhi) = diagrama(numaa+1, diagPsi, diagPhi) +1

			totalpoints(aa) = totalpoints(aa) +1
			totalpoints(numaa+1) = totalpoints(numaa+1) +1
		end do
	end do
	close(10)
	print '("Analizado: [", i3, "/", i3, "]")', fileind, numfiles
end do

do waa=1, numaa+1
	do diagrow=1, resol
		do diagcol=1, resol
			if (diagrama(waa, diagrow, diagcol).eq.0) then
				diagrama(waa, diagrow, diagcol) = -10
			else
				diagrama(waa, diagrow, diagcol) = (diagrama(waa, diagrow, diagcol)/totalpoints(waa))*100
			end if
		end do
	end do
end do
	
! Guardado de resultados

! Global
do diagrow=1, resol
	write(20,*) (diagrama(numaa+1, diagrow, diagcol), diagcol=1, resol)		
end do

! Individual
do waa=1, numaa
	do diagrow=1, resol
		write(outputs(waa),*) (diagrama(waa, diagrow, diagcol), diagcol=1, resol)	
	end do
end do

end program ramachandran




