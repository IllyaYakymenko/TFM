    
!	#################################################################################
!	#
!	#   Obtención del mapa de densidades de los desplazamientos
!	#   
!	#   Input: Channel 10 (fichero _raman.dat)
!	#   Output: Channel 20: Mapa de densidades DeltaPhi/DeltaPsi
!	#	          Channel 30: Mapa de densidades DeltaPhi+DeltaPsi/DeltaPhi-DeltaPsi
!	#
!	#################################################################################

program dens_desplazamientos

implicit none

! Variables de input

integer numfiles, nrow, ncol, PhiRel, PsiRel
real PhiCenter, PsiCenter, PhiRango, PsiRango, resol
character*100, allocatable :: filenames(:)

! Variables de análisis

double precision, allocatable :: diedros(:), diedrosant(:), desps(:)
double precision dPhi, dPsi, infPhi, supPhi, infPsi, supPsi, dsuma, dresta, N, Phi, Psi
integer allRango, limRango, numaa, to, from, regaa
double precision, allocatable :: densPhiPsi(:,:), denssumaresta(:,:)
integer posdPhi, posdPsi, posdsuma, posdresta

! Variables de iteración

integer fileind, r, aa, ind1, ind2, rowiter

read(5,*) numfiles
allocate(filenames(numfiles))
read(5,*) filenames
read(5,*) ncol
read(5,*) nrow
read(5,*) PhiCenter, PsiCenter, PhiRango, PsiRango
read(5,*) resol
read(5,*) limRango
read(5,*) PhiRel, PsiRel

allRango = int((limRango*2) / resol)
numaa = nrow/2

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

allocate(diedros(nrow))
allocate(diedrosant(nrow))
allocate(desps(nrow))
allocate(densPhiPsi(allRango, allRango))
allocate(denssumaresta(allRango, allRango))

! Delimitación de la región a analizar

infPhi = (PhiCenter - PhiRango)
supPhi = (PhiCenter + PhiRango)
infPsi = (PsiCenter - PhiRango)
supPsi = (PsiCenter + PsiRango)

if(infPhi.lt.-180) infPhi = infPhi+360
if(supPhi.gt.180) supPhi = supPhi-360
if(infPsi.lt.-180) infPsi = infPsi+360
if(supPsi.gt.180) supPsi = supPsi-360

print '("Registrando datos entre Phi(", f7.2, ":", f7.2, "); Psi(", f7.2, ":", f7.2, ")")', infPhi, supPhi, infPsi, supPsi

! Matrices vacias (nulas)

do ind1=1, allRango
	do ind2=1, allRango
		densPhiPsi(ind1, ind2) = 0.d0
		denssumaresta(ind1, ind2) = 0.d0
	end do
end do

N = 0

print*, "Registrando valores:"

do fileind=1, numfiles
	open(10, file = filenames(fileind))
	read(10,*) diedrosant

	do rowiter=2, ncol
		read(10,*) diedros
		
		do r=1, nrow
			desps(r) = (diedros(r) - diedrosant(r)) - 1*(int((diedros(r) - diedrosant(r))/180))*360
		end do

		do aa=from, to
			Phi = diedros(2*aa-1 + 2*PhiRel)
			Psi = diedros(2*aa + 2*PsiRel)

			! Si el aminoácido se encuentra en el límite de conformaciones especificado se efectúan los cálculos
			! Se registra la densidad como sumatorio para la posición correspondiente
			if((Phi.ge.infPhi).and.(Phi.le.supPhi).and.(Psi.ge.infPsi).and.(Psi.le.supPsi)) then
				N = N+1

				dPhi = desps(2*aa-1 + 2*PhiRel)
				dPsi = desps(2*aa + 2*PsiRel)
				dsuma = dPhi+dPsi
				dresta = dPhi-dPsi

				posdPhi = int((dPhi + limRango)/resol)
				posdPsi = int((dPsi + limRango)/resol)
				posdsuma = int((dsuma + limRango)/resol)
				posdresta = int((dresta + limRango)/resol)

				densPhiPsi(posdPsi, posdPhi) = densPhiPsi(posdPsi, posdPhi) +1
				denssumaresta(posdsuma, posdresta) = denssumaresta(posdsuma, posdresta) +1

			end if
		
		end do
		diedrosant=diedros
	
	end do
	print '("Registrando valores [", i3,"/", i3, "]")', fileind, numfiles
	close(10)
end do

! Normalización (sobre 100) de las densidades

do ind1=1, allRango
	do ind2=1, allRango
		densPhiPsi(ind1, ind2) = (densPhiPsi(ind1, ind2)/N)*100
		denssumaresta(ind1, ind2) = (denssumaresta(ind1, ind2)/N)*100
	end do
end do

! Escritura de los mapas

do ind1=1, allRango
	write(20, *) (densPhiPsi(ind1, ind2), ind2=1, allRango) 
end do

do ind1=1, allRango
	write(30, *) (denssumaresta(ind1, ind2), ind2=1, allRango) 
end do


end program dens_desplazamientos
