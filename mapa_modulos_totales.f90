! ############################################################################
! #
! #   Obtención mapas de módulos para cada aminoácido
! #   
! #   Input: Channel 10 (fichero _raman.dat)
! #   Output: Channel 20 -> 20 + número de aminoácidos:
! #                  Mapas del valor medio de los desplazamientos de T
! #           Channel 50 -> 50 + número de aminoácidos:
! #                  Mapas del valor medio de los desplazamientos de T^2
! #
! #
! ############################################################################

program mapa_modulos_totales

implicit none

! Datos de entrada (archivo.in)

integer :: numfiles, resol, nrow, ncol, PhiRel, PsiRel
character*200, allocatable :: filenames(:)
real :: PhiMin, PhiMax, PsiMin, PsiMax

! Definición de la matriz del diagrama

double precision, allocatable :: sumaT(:,:,:), sumaTCuad(:,:,:), mapaT(:,:,:), mapaTCuad(:,:,:)
double precision, allocatable :: N(:,:,:)
double precision :: rangoPhi, rangoPsi, dPhi, dPsi, modT, modTCuad
integer :: diagPhi, diagPsi, numaa, from, to
double precision, allocatable :: diedros(:), diedrosant(:), desps(:)
integer, allocatable :: outputsT(:), outputsTCuad(:)

! Variables de control

integer ch, row, col, fileind, rowiter, dr, dc, aa, i, diagaa, regaa, di, waa

read(5,*) numfiles

allocate(filenames(numfiles))

read(5,*) filenames
read(5,*) nrow
read(5,*) ncol
read(5,*) PhiMin, PhiMax, PsiMin, PsiMax
read(5,*) resol
read(5,*) PhiRel, PsiRel

numaa = ncol/2

! Evaluación de los aminoácidos a analizar
! Su valor depende de la relatividad de los ángulos
! En casos de PhiRel y PsiRel distinta a 0 hay aminoácidos que no se analizarían

from = 1
to = numaa

if((PhiRel.lt.0).or.(PsiRel.lt.0)) then
	from = 1- min(PhiRel, PsiRel)
end if

if((PhiRel.gt.0).or.(PsiRel.gt.0)) then
	to = numaa - max(PhiRel, PsiRel)
end if

regaa = to - from +1

allocate(sumaT(from:to+1, resol, resol))
allocate(sumaTCuad(from:to+1, resol, resol))
allocate(N(from:to+1,resol,resol))
allocate(mapaT(from:to+1,resol,resol))
allocate(mapaTCuad(from:to+1,resol,resol))
allocate(diedros(ncol))
allocate(diedrosant(ncol))
allocate(desps(ncol))
allocate(outputsT(from:to))
allocate(outputsTCuad(from:to))

! Creación de canales de output
! En outputs se guardan por orden las salidas de cada amonoácido
! Hay tantos outputs como aminoácidos en el polipéptido
! Las frecuencias de todos los aminoácidos se guardan en la posición
! numaa+1 de la primera dimensión de diagrama

do ch=from, to
	outputsT(ch) = 20+ch
	outputsTCuad(ch) = 50+ch
end do

! Matrices vacias (nulas)
do diagaa=from, to+1
	do row=1, resol
		do col=1, resol
			sumaT(diagaa, row, col) = 0.d0
			sumaTCuad(diagaa, row, col) = 0.d0
			N(diagaa,row, col) = 0.d0
			mapaT(diagaa, row, col) = 0.d0
			mapaTCuad(diagaa, row, col) = 0.d0
		end do
	end do
end do

rangoPhi = (PhiMax - PhiMin)/resol
rangoPsi = (PsiMax - PsiMin)/resol

print*, 'Calculando mapas de módulos totales'

do fileind=1, numfiles
	open(10, file = filenames(fileind))
	read(10,*) diedrosant

	do rowiter=2, nrow
		read(10,*) diedros

		do di=1, ncol
		desps(di) = (diedros(di) - diedrosant(di))-1*int((diedros(di) - diedrosant(di))/180)*360	
		end do

		do aa=from, to
			! Posición en las matrices
			diagPhi = int((diedros(2*aa-1)+PhiMax)/rangoPhi)+1
			diagPsi = int((diedros(2*aa)+PsiMax)/rangoPsi)+1
			
			if(diagPhi.ge.resol) then
				diagPhi=resol
			end if

			if(diagPsi.ge.resol) then
				diagPsi=resol
			end if

			! Desplazamiento del ángulo solicitado
			dPhi = desps(2*aa-1 + 2*PhiRel)
			dPsi = desps(2*aa + 2*PsiRel)

			! Cálculo de los módulos de interés
			modT = sqrt(dPhi*dPhi + dPsi*dPsi)
			modTCuad = dPhi*dPhi + dPsi*dPsi

			sumaT(aa, diagPsi, diagPhi) = sumaT(aa, diagPsi, diagPhi) + modT
			sumaTCuad(aa, diagPsi, diagPhi) = sumaTCuad(aa, diagPsi, diagPhi) + modTCuad
			N(aa, diagPsi, diagPhi) = N(aa, diagPsi, diagPhi) + 1

			sumaT(to+1, diagPsi, diagPhi) = sumaT(to+1, diagPsi, diagPhi) + modT
			sumaTCuad(to+1, diagPsi, diagPhi) = sumaTCuad(to+1, diagPsi, diagPhi) + modTCuad
			N(to+1, diagPsi, diagPhi) = N(to+1, diagPsi, diagPhi) + 1
		
		end do
		diedrosant = diedros
	end do
	close(10)
	print '("Analizado: [", i3, "/", i3, "]")', fileind, numfiles
end do

! Cálculo de los valores medios de los desplazamientos en cada posición del mapa
do aa=from, to+1
	do dr=1, resol
		do dc=1, resol
			if (N(aa, dr, dc).gt.0) then
				mapaT(aa, dr, dc) = sumaT(aa, dr, dc)/N(aa, dr, dc)
				mapaTCuad(aa, dr, dc) = sumaTCuad(aa, dr, dc)/N(aa, dr, dc)
			else
				mapaT(aa, dr, dc) = -10000
				mapaTCuad(aa, dr, dc) = -10000
			end if
		end do
	end do
end do
	
! Guardado de resultados

! Global
do dr=1, resol
	write(20,*) (mapaT(to+1, dr, dc), dc=1, resol)		
end do

! Individual
do waa=from, to
	do dr=1, resol
		write(outputsT(waa),*) (mapaT(waa, dr, dc), dc=1, resol)	
	end do
end do


! Global
do dr=1, resol
	write(50,*) (mapaTCuad(to+1, dr, dc), dc=1, resol)		
end do

! Individual
do waa=from, to
	do dr=1, resol
		write(outputsTCuad(waa),*) (mapaTCuad(waa, dr, dc), dc=1, resol)	
	end do
end do

end program mapa_modulos_totales




