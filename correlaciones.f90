!	############################################################################
!	#
!	#   Obtención del diagrama de Ramachandran para cada
!	#   aminoácido.
!	#   Input: Channel 10 (fichero _raman.dat)
!	#   Output: Channel 20: Mapa de correlaciones de todos los aminoácidos
!	#	        21-> 20 + número de aminoácidos:
!	#	        Mapa de correlaciones de cada aminoácido
!	#              Canales dependientes de los aminoácidos de análisis
!	#
!	############################################################################


program mapa_correlaciones

implicit none

! Datos de entrada (archivo.in)

integer :: numfiles, resol, numrow, numcol, PhiRelAA, PsiRelAA
character*200, allocatable :: filenames(:)
real :: PhiMin, PhiMax, PsiMin, PsiMax

! Definición de las matrices de análisis

double precision, allocatable :: correlaciones(:,:,:), SPhiPsi(:,:,:), SPhi(:,:,:) 
double precision, allocatable :: SPsi(:,:,:), SPhiCuad(:,:,:), SPsiCuad(:,:,:) 
integer, allocatable :: N(:,:,:)
double precision :: rangoPhi, rangoPsi, despPhi, despPsi, dSPhi, dSPsi, dSPhiPsi
integer :: diagPhi, diagPsi, numaa
double precision, allocatable :: diedros(:), diedrosant(:), desplazamientos(:)
integer, allocatable :: outputs(:)

! Variables de control

integer ch, row, col, fileind, rowiter, dr, dc, aa, i, diagaa, waa, from, to, regaa, di

read(5,*) numfiles

allocate(filenames(numfiles))

read(5,*) filenames
read(5,*) numrow
read(5,*) numcol
read(5,*) PhiMin, PhiMax, PsiMin, PsiMax
read(5,*) PhiRelAA, PsiRelAA
read(5,*) resol

numaa = numcol/2

! Evaluación de los aminoácidos a analizar
! Su valor depende de la relatividad de los ángulos
! En casos de PhiRel y PsiRel distinta a 0 hay aminoácidos que no se analizarían
from = 1
to = numaa

if ((PhiRelAA.lt.0).or.(PsiRelAA.lt.0)) then
	from = 1 - min(PhiRelAA, PsiRelAA)
end if

if ((PhiRelAA.gt.0).or.(PsiRelAA.gt.0)) then
	to = numaa - max(PhiRelAA, PsiRelAA)
end if

regaa = to - from +1

allocate(correlaciones(from:to+1,resol,resol))
allocate(SPhiPsi(from:to+1,resol,resol))
allocate(SPhi(from:to+1,resol,resol))
allocate(SPsi(from:to+1,resol,resol))
allocate(SPhiCuad(from:to+1,resol,resol))
allocate(SPsiCuad(from:to+1,resol,resol))
allocate(N(from:to+1,resol,resol))
allocate(diedros(numcol))
allocate(diedrosant(numcol))
allocate(desplazamientos(numcol))
allocate(outputs(from:to))

! Creación de canales de output
! En outputs se guardan por orden las salidas de cada amonoácido
! Hay tantos outputs como aminoácidos en el polipéptido
! Las frecuencias de todos los aminoácidos se guardan en la posición

do ch=from, to
	outputs(ch) = 20+ch	
end do

! Matrices vacias (nulas)
do diagaa=from, to+1
	do row=1, resol
		do col=1, resol
			correlaciones(diagaa,row, col) = 0.d0
			SPhiPsi(diagaa,row, col) = 0.d0
			SPhi(diagaa,row, col) = 0.d0
			SPsi(diagaa,row, col) = 0.d0
			SPhiCuad(diagaa,row, col) = 0.d0
			SPsiCuad(diagaa,row, col) = 0.d0
			N(diagaa,row, col) = 0.d0
		end do
	end do
end do

rangoPhi = (PhiMax - PhiMin)/resol
rangoPsi = (PsiMax - PsiMin)/resol

print*, 'Calculando mapas de correlaciones'

do fileind=1, numfiles
	open(10, file = filenames(fileind))
	read(10,*) diedrosant

	do rowiter=2, numrow
		read(10,*) diedros

		do di=1, numcol
		desplazamientos(di) = (diedros(di) - diedrosant(di))-1*int((diedros(di) - diedrosant(di))/180)*360	
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
			despPhi = desplazamientos(2*aa+2*PhiRelAA-1)
			despPsi = desplazamientos(2*aa+2*PsiRelAA)

			! Registro de los desplazamientos en las distintas matrices
			SPhiPsi(aa, diagPsi, diagPhi) = SPhiPsi(aa, diagPsi, diagPhi) + despPhi*despPsi
			SPhi(aa, diagPsi, diagPhi) = SPhi(aa, diagPsi, diagPhi) + despPhi
			SPsi(aa, diagPsi, diagPhi) = SPsi(aa, diagPsi, diagPhi) + despPsi
			SPhiCuad(aa, diagPsi, diagPhi) = SPhiCuad(aa, diagPsi, diagPhi) + despPhi*despPhi
			SPsiCuad(aa, diagPsi, diagPhi) = SPsiCuad(aa, diagPsi, diagPhi) + despPsi*despPsi
			N(aa, diagPsi, diagPhi) = N(aa, diagPsi, diagPhi) + 1

			SPhiPsi(to+1, diagPsi, diagPhi) = SPhiPsi(to+1, diagPsi, diagPhi) + despPhi*despPsi
			SPhi(to+1, diagPsi, diagPhi) = SPhi(to+1, diagPsi, diagPhi) + despPhi
			SPsi(to+1, diagPsi, diagPhi) = SPsi(to+1, diagPsi, diagPhi) + despPsi
			SPhiCuad(to+1, diagPsi, diagPhi) = SPhiCuad(to+1, diagPsi, diagPhi) + despPhi*despPhi
			SPsiCuad(to+1, diagPsi, diagPhi) = SPsiCuad(to+1, diagPsi, diagPhi) + despPsi*despPsi
			N(to+1, diagPsi, diagPhi) = N(to+1, diagPsi, diagPhi) + 1
			
		end do
		diedrosant = diedros
	end do
	close(10)
	print '("Analizado: [", i3, "/", i3, "]")', fileind, numfiles
end do

! Cálculo de correlaciones

do waa=from, to+1
	do dr=1, resol
		do dc=1, resol
			if (N(waa, dr, dc).gt.1) then
				dSPhiPsi = dfloat(N(waa, dr, dc))*SPhiPsi(waa, dr, dc) - SPhi(waa, dr, dc)*SPsi(waa, dr, dc)
				dSPhi = sqrt(dfloat(N(waa, dr, dc))*SPhiCuad(waa, dr, dc) - SPhi(waa, dr, dc)*SPhi(waa, dr, dc))
				dSPsi = sqrt(dfloat(N(waa, dr, dc))*SPsiCuad(waa, dr, dc) - SPsi(waa, dr, dc)*SPsi(waa, dr, dc))

				correlaciones(waa, dr, dc) = dSPhiPsi / (dSPhi*dSPsi)
			else
			! Si N es 0 o 1, no existe correlación, por lo cual, esos casos no aparecen en el mapa
				correlaciones(waa, dr, dc) = -10
			end if
		end do
	end do
end do
	
! Guardado de resultados

! Global
do dr=1, resol
	write(20,*) (correlaciones(to+1, dr, dc), dc=1, resol)		
end do

! Individual
do waa=from, to
	do dr=1, resol
		write(outputs(waa),*) (correlaciones(waa, dr, dc), dc=1, resol)	
	end do
end do

end program mapa_correlaciones




