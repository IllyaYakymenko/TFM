program distribuciones

implicit none

! Variables de input

integer numfiles, nrow, ncol, PhiRel, PsiRel
real PhiCenter, PsiCenter, PhiRango, PsiRango, resol
character*100, allocatable :: filenames(:)

! Variables de análisis

double precision, allocatable :: diedros(:), diedrosant(:), desps(:)
integer numaa, posPhi, posPsi, posT, posTs, posTCuad, possuma, posresta, posTCuadsigno
real infPhi, infPsi, supPhi, supPsi
integer, allocatable :: rangodist(:)
double precision, allocatable :: limites(:,:)
double precision dPhi, dPsi, dT, dTs, dTCuadsigno, dTCuad, dsuma, dresta, Phi, Psi, liminf, limsup
double precision, allocatable :: distPhi(:,:,:), distPsi(:,:,:), distT(:,:,:), distTs(:,:,:)
double precision, allocatable :: distTCuad(:,:,:), distsuma(:,:,:), distresta(:,:,:)
double precision, allocatable :: distTCuadSigno(:,:,:)
double precision, allocatable :: sumatorio(:,:), medias(:,:), sumaG(:,:), G(:,:), N(:)
double precision, allocatable :: sumaTsigno(:,:), mediaTsigno(:,:), sumaGsigno(:,:), Gsigno(:,:), Nsigno(:,:)
double precision, allocatable :: sumaTCuadsigno(:,:), mediaTCuadsigno(:,:), sumaGCuadsigno(:,:), GCuadsigno(:,:)
character*6 glob

! Variables de iteración

integer fileind, r, aa, ind1, ind2, rowiter, sclas, from, to, regaa

read(5,*) numfiles
allocate(filenames(numfiles))
read(5,*) filenames
read(5,*) ncol
read(5,*) nrow
read(5,*) PhiCenter, PsiCenter, PhiRango, PsiRango
read(5,*) resol
read(5,*) PhiRel, PsiRel

! resolTCuad = resol*10 ! Valor del límite muy elevado para el histograma

numaa=nrow/2

from = 1
to = numaa

if((PhiRel.lt.0).or.(PsiRel.lt.0)) then
	from = 1- min(PhiRel, PsiRel)
end if

if((PhiRel.gt.0).or.(PsiRel.gt.0)) then
	to = numaa - max(PhiRel, PsiRel)
end if

regaa = to - from +1

allocate(limites(8,2))
allocate(rangodist(8))
allocate(sumatorio(from:to+1,8))
allocate(medias(from:to+1,8))
allocate(sumaG(from:to+1,8))
allocate(G(from:to+1,8))
allocate(N(from:to+1))
allocate(sumaTsigno(from:to+1,2))
allocate(mediaTsigno(from:to+1,2))
allocate(sumaGsigno(from:to+1,2))
allocate(Gsigno(from:to+1,2))
allocate(sumaTCuadsigno(from:to+1,2))
allocate(mediaTCuadsigno(from:to+1,2))
allocate(sumaGCuadsigno(from:to+1,2))
allocate(GCuadsigno(from:to+1,2))
allocate(Nsigno(from:to+1,2))
allocate(diedros(nrow))
allocate(diedrosant(nrow))
allocate(desps(nrow))

! Cuando se establece un límite de 8 en una dimensión, ese valor corresponde a:
! 1 = Phi
! 2 = Psi
! 3 = T (módulo)
! 4 = Ts (módulo dependiente de signos)
! 5 = TCuad (módulo^2)
! 6 = Phi+Phi (suma)
! 7 = Phi-Psi (resta)
! 8 = TCuadSigno

infPhi = (PhiCenter - PhiRango)
supPhi = (PhiCenter + PhiRango)
infPsi = (PsiCenter - PhiRango)
supPsi = (PsiCenter + PsiRango)

if(infPhi.lt.-180) infPhi = infPhi+360
if(supPhi.gt.180) supPhi = supPhi-360
if(infPsi.lt.-180) infPsi = infPsi+360
if(supPsi.gt.180) supPsi = supPsi-360

! print*, infPhi, supPhi, infPsi, supPsi

print '("Registrando datos entre Phi(", f7.2, ":", f7.2, "); Psi(", f7.2, ":", f7.2, ")")', infPhi, supPhi, infPsi, supPsi

do ind1=1, 8
	do ind2=1, 2
		limites(ind1, ind2) = 0.d0
		
		do aa=from, to+1
			sumatorio(aa,ind1) = 0.d0
			sumaG(aa,ind1) = 0.d0
		end do
	end do
end do

do aa=from, to+1
	N(aa) = 0.d0

	do ind1=1, 2
		Nsigno(aa, ind1) = 0.d0
		sumaTsigno(aa,ind1) = 0.d0
		sumaGsigno(aa,ind1) = 0.d0
		sumaTCuadsigno(aa,ind1) = 0.d0
		sumaGCuadsigno(aa,ind1) = 0.d0
	end do
end do

! Búsqueda de límites y medias

print*, "Buscando límites y calculando medias:"

do fileind=1, numfiles
	open(10, file = filenames(fileind))
	read(10,*) diedrosant

	do rowiter=2, ncol
		read(10,*) diedros
		
		do r=1, nrow
			desps(r) = (diedros(r) - diedrosant(r)) - 1*(int((diedros(r) - diedrosant(r))/180))*360
		end do

		do aa=from, to
			Phi = diedros((2*aa-1)+2*PhiRel)
			Psi = diedros(2*aa + 2*PsiRel)

			if((Phi.ge.infPhi).and.(Phi.le.supPhi).and.(Psi.ge.infPsi).and.(Psi.le.supPsi)) then
				N(aa) = N(aa)+1
				N(to+1) = N(to+1)+1

				dPhi = desps((2*aa-1) +2*PhiRel)
				dPsi = desps(2*aa + 2*PsiRel)
				dT = sqrt(dPhi*dPhi + dPsi*dPsi)
				dTCuad = dPhi*dPhi + dPsi*dPsi
				dsuma = dPhi + dPsi
				dresta = dPhi - dPsi
				dTs = dT
				dTCuadsigno = dTCuad

				sclas = 1
				if(((dPhi.gt. 0.d0).and.(dPsi.lt. 0.d0)).or.((dPhi.lt. 0.d0).and.(dPsi.gt. 0.d0))) then
					dTs = dTs*(-1)
					dTCuadsigno = dTCuadsigno*(-1)
					sclas = 2
				end if

				Nsigno(aa, sclas) = Nsigno(aa, sclas) +1
				Nsigno(to+1, sclas) = Nsigno(to+1, sclas) +1

				limites(1,1) = min(limites(1,1), dPhi)
				limites(1,2) = max(limites(1,2), dPhi)
				sumatorio(aa,1) = sumatorio(aa,1) + dPhi
				sumatorio(to+1,1) = sumatorio(to+1,1) + dPhi

				limites(2,1) = min(limites(2,1), dPsi)
				limites(2,2) = max(limites(2,2), dPsi)
				sumatorio(aa,2) = sumatorio(aa,2) + dPsi
				sumatorio(to+1,2) = sumatorio(to+1,2) + dPsi

				limites(3,1) = min(limites(3,1), dT)
				limites(3,2) = max(limites(3,2), dT)
				sumatorio(aa,3) = sumatorio(aa,3) + dT
				sumatorio(to+1,3) = sumatorio(to+1,3) + dT

				limites(4,1) = min(limites(4,1), dTs)
				limites(4,2) = max(limites(4,2), dTs)
				sumatorio(aa,4) = sumatorio(aa,4) + dTs
				sumatorio(to+1,4) = sumatorio(to+1,4) + dTs

				limites(5,1) = min(limites(5,1), dTCuad)
				limites(5,2) = max(limites(5,2), dTCuad)
				sumatorio(aa,5) = sumatorio(aa,5) + dTCuad
				sumatorio(to+1,5) = sumatorio(to+1,5) + dTCuad

				limites(6,1) = min(limites(6,1), dsuma)
				limites(6,2) = max(limites(6,2), dsuma)
				sumatorio(aa,6) = sumatorio(aa,6) + dsuma
				sumatorio(to+1,6) = sumatorio(to+1,6) + dsuma

				limites(7,1) = min(limites(7,1), dresta)
				limites(7,2) = max(limites(7,2), dresta)
				sumatorio(aa,7) = sumatorio(aa,7) + dresta
				sumatorio(to+1,7) = sumatorio(to+1,7) + dresta

				limites(8,1) = min(limites(8,1), dTCuadsigno)
				limites(8,2) = max(limites(8,2), dTCuadsigno)
				sumatorio(aa,8) = sumatorio(aa,8) + dTCuadsigno
				sumatorio(to+1,8) = sumatorio(to+1,8) + dTCuadsigno

				sumaTsigno(aa, sclas) = sumaTsigno(aa, sclas) + dT
				sumaTsigno(to+1, sclas) = sumaTsigno(to+1, sclas) + dT
				
				sumaTCuadsigno(aa, sclas) = sumaTCuadsigno(aa, sclas) + dTCuad
				sumaTCuadsigno(to+1, sclas) = sumaTCuadsigno(to+1, sclas) + dTCuad

			end if
		
		end do
		diedrosant=diedros
	
	end do
	print '("Buscando límites y calculando medias [", i3,"/", i3, "]")', fileind, numfiles
	close(10)
end do

! Establecimiento de los límites de cada histograma

do ind1=1, 8
	liminf = limites(ind1, 1)
	limsup = limites(ind1, 2)

	do aa=from, to+1
		medias(aa,ind1) = sumatorio(aa,ind1)/N(aa)
	end do

	rangodist(ind1) = int(((int(limsup)+1)-(int(liminf)-1))/resol) + 1
end do
	
! rangodist(5) = int(((int(limites(5,2))+1)-(int(limites(5,1))-1))/resol) + 1

do aa=from, to+1
	do ind1=1, 2
		mediaTsigno(aa, ind1) = sumaTsigno(aa,ind1)/Nsigno(aa, ind1)
		mediaTCuadsigno(aa, ind1) = sumaTCuadsigno(aa,ind1)/Nsigno(aa, ind1)
	end do
end do

allocate(distPhi(from:to+1, 2, rangodist(1)))
allocate(distPsi(from:to+1, 2, rangodist(2)))
allocate(distT(from:to+1, 2, rangodist(3)))
allocate(distTs(from:to+1, 2, rangodist(4)))
allocate(distTCuad(from:to+1, 2, rangodist(5)))
allocate(distsuma(from:to+1, 2, rangodist(6)))
allocate(distresta(from:to+1, 2, rangodist(7)))
allocate(distTCuadsigno(from:to+1, 2, rangodist(8)))

! Iniciación de los histogramas de distribución

do aa=from, to+1
	do ind1=1, rangodist(1)
		liminf = limites(1,1)
		distPhi(aa, 1, ind1) = (int(liminf)-1) + (ind1-1)*resol
		distPhi(aa, 2, ind1) = 0.d0
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(2)
		liminf = limites(2,1)
		distPsi(aa, 1, ind1) = (int(liminf)-1) + (ind1-1)*resol
		distPsi(aa, 2, ind1) = 0.d0
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(3)
		liminf = limites(3,1)
		distT(aa, 1, ind1) = (int(liminf)-1) + (ind1-1)*resol
		distT(aa, 2, ind1) = 0.d0
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(4)
		liminf = limites(4,1)
		distTs(aa, 1, ind1) = (int(liminf)-1) + (ind1-1)*resol
		distTs(aa, 2, ind1) = 0.d0
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(5)
		liminf = limites(5,1)
		distTCuad(aa, 1, ind1) = (int(liminf)-1) + (ind1-1)*resol
		distTCuad(aa, 2, ind1) = 0.d0
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(6)
		liminf = limites(6,1)
		distsuma(aa, 1, ind1) = (int(liminf)-1) + (ind1-1)*resol
		distsuma(aa, 2, ind1) = 0.d0
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(7)
		liminf = limites(7,1)
		distresta(aa, 1, ind1) = (int(liminf)-1) + (ind1-1)*resol
		distresta(aa, 2, ind1) = 0.d0
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(8)
		liminf = limites(8,1)
		distTCuadsigno(aa, 1, ind1) = (int(liminf)-1) + (ind1-1)*resol
		distTCuadsigno(aa, 2, ind1) = 0.d0
	end do
end do

! Registro de histogramas y cálculo de desviación estándar

print*, "Registrando histogramas y calculando desviación estándar:"

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

			if((Phi.ge.infPhi).and.(Phi.le.supPhi).and.(Psi.ge.infPsi).and.(Psi.le.supPsi)) then
				dPhi = desps(2*aa-1 + 2*PhiRel)
				dPsi = desps(2*aa + 2*PsiRel)
				dT = sqrt(dPhi*dPhi + dPsi*dPsi)
				dTCuad = dPhi*dPhi + dPsi*dPsi
				dsuma = dPhi + dPsi
				dresta = dPhi - dPsi
				dTs = dT
				dTCuadsigno = dTCuad

				sclas = 1
				if(((dPhi.gt. 0.d0).and.(dPsi.lt. 0.d0)).or.((dPhi.lt. 0.d0).and.(dPsi.gt. 0.d0))) then
					dTs = dTs*(-1)
					dTCuadsigno = dTCuadsigno*(-1)
					sclas = 2
				end if
				
				posPhi = int((dPhi - (int(limites(1,1))-1))/resol) +1
				distPhi(aa, 2, posPhi) = distPhi(aa, 2, posPhi) +1
				distPhi(to+1, 2, posPhi) = distPhi(to+1, 2, posPhi) +1
				sumaG(aa,1) = sumaG(aa,1) + (medias(aa,1)-dPhi)**2
				sumaG(to+1,1) = sumaG(to+1,1) + (medias(to+1,1)-dPhi)**2

				posPsi = int((dPsi - (int(limites(2,1))-1))/resol) +1
				distPsi(aa, 2, posPsi) = distPsi(aa, 2, posPsi) +1
				distPsi(to+1, 2, posPsi) = distPsi(to+1, 2, posPsi) +1
				sumaG(aa,2) = sumaG(aa,2) + (medias(aa,2)-dPsi)**2
				sumaG(to+1,2) = sumaG(to+1,2) + (medias(to+1,2)-dPsi)**2

				posT = int((dT - (int(limites(3,1))-1))/resol) +1
				distT(aa, 2, posT) = distT(aa, 2, posT) +1
				distT(to+1, 2, posT) = distT(to+1, 2, posT) +1
				sumaG(aa,3) = sumaG(aa,3) + (medias(aa,3)-dT)**2
				sumaG(to+1,3) = sumaG(to+1,3) + (medias(to+1,3)-dT)**2

				posTs = int((dTs - (int(limites(4,1))-1))/resol) +1
				distTs(aa, 2, posTs) = distTs(aa, 2, posTs) +1
				distTs(to+1, 2, posTs) = distTs(to+1, 2, posTs) +1
				sumaG(aa,4) = sumaG(aa,4) + (medias(aa,4)-dTs)**2
				sumaG(to+1,4) = sumaG(to+1,4) + (medias(to+1,4)-dTs)**2

				posTCuad = int((dTCuad - (int(limites(5,1))-1))/resol) +1
				distTCuad(aa, 2, posTCuad) = distTCuad(aa, 2, posTCuad) +1
				distTCuad(to+1, 2, posTCuad) = distTCuad(to+1, 2, posTCuad) +1
				sumaG(aa,5) = sumaG(aa,5) + (medias(aa,5)-dTCuad)**2
				sumaG(to+1,5) = sumaG(to+1,5) + (medias(to+1,5)-dTCuad)**2

				possuma = int((dsuma - (int(limites(6,1))-1))/resol) +1
				distsuma(aa, 2, possuma) = distsuma(aa, 2, possuma) +1
				distsuma(to+1, 2, possuma) = distsuma(to+1, 2, possuma) +1
				sumaG(aa,6) = sumaG(aa,6) + (medias(aa,6)-dsuma)**2
				sumaG(to+1,6) = sumaG(to+1,6) + (medias(to+1,6)-dsuma)**2

				posresta = int((dresta - (int(limites(7,1))-1))/resol) +1
				distresta(aa, 2, posresta) = distresta(aa, 2, posresta) +1
				distresta(to+1, 2, posresta) = distresta(to+1, 2, posresta) +1
				sumaG(aa,7) = sumaG(aa,7) + (medias(aa,7)-dresta)**2
				sumaG(to+1,7) = sumaG(to+1,7) + (medias(to+1,7)-dresta)**2

				posTCuadsigno = int((dTCuadsigno - (int(limites(8,1))-1))/resol) +1
				distTCuadsigno(aa, 2, posTCuadsigno) = distTCuadsigno(aa, 2, posTCuadsigno) +1
				distTCuadsigno(to+1, 2, posTCuadsigno) = distTCuadsigno(to+1, 2, posTCuadsigno) +1
				sumaG(aa,8) = sumaG(aa,8) + (medias(aa,8)-dTCuadsigno)**2
				sumaG(to+1,8) = sumaG(to+1,8) + (medias(to+1,8)-dTCuadsigno)**2

				sumaGsigno(aa, sclas) = sumaGsigno(aa, sclas) + (mediaTsigno(aa, sclas)-dT)**2
				sumaGsigno(to+1, sclas) = sumaGsigno(to+1, sclas) + (mediaTsigno(to+1, sclas)-dT)**2

				sumaGCuadsigno(aa, sclas) = sumaGCuadsigno(aa, sclas) + (mediaTCuadsigno(aa, sclas)-dTCuad)**2
				sumaGCuadsigno(to+1, sclas) = sumaGCuadsigno(to+1, sclas) + (mediaTCuadsigno(to+1, sclas)-dTCuad)**2

			end if
		
		end do
		diedrosant=diedros
	
	end do
	print '("Registrando histogramas y calculando desviación estándar [", i3,"/", i3, "]")', fileind, numfiles
	close(10)
end do

do aa=from, to+1
	do ind1=1, 8
		G(aa,ind1) = sqrt(sumaG(aa,ind1)/N(aa))
	end do
end do

do aa=from, to+1
	do ind1=1, 2
		Gsigno(aa, ind1) = sqrt(sumaGsigno(aa, ind1)/Nsigno(aa, ind1))
		GCuadsigno(aa, ind1) = sqrt(sumaGCuadsigno(aa, ind1)/Nsigno(aa, ind1))
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(1)
		distPhi(aa, 2, ind1) = (distPhi(aa, 2, ind1)*100)/N(aa)
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(2)
		distPsi(aa, 2, ind1) = (distPsi(aa, 2, ind1)*100)/N(aa)
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(3)
		distT(aa, 2, ind1) = (distT(aa, 2, ind1)*100)/N(aa)
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(4)
		distTs(aa, 2, ind1) = (distTs(aa, 2, ind1)*100)/N(aa)
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(5)
		distTCuad(aa, 2, ind1) = (distTCuad(aa, 2, ind1)*100)/N(aa)
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(6)
		distsuma(aa, 2, ind1) = (distsuma(aa, 2, ind1)*100)/N(aa)
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(7)
		distresta(aa, 2, ind1) = (distresta(aa, 2, ind1)*100)/N(aa)
	end do
end do

do aa=from, to+1
	do ind1=1, rangodist(8)
		distTCuadsigno(aa, 2, ind1) = (distTCuadsigno(aa, 2, ind1)*100)/N(aa)
	end do
end do


! Escribiendo los datos en ficheros

glob = "Global"

!Phi
do aa=from, to
	do ind1=1, rangodist(1)
		write(21,"(i6, 3x, f13.8, 3x, f13.8)") aa, distPhi(aa, 1, ind1), distPhi(aa, 2, ind1)
	end do
end do

do ind1=1, rangodist(1)
	write(21,"(a6, 3x, f13.8, 3x, f13.8)") glob, distPhi(to+1, 1, ind1), distPhi(to+1, 2, ind1)
end do

!Psi
do aa=from, to
	do ind1=1, rangodist(2)
		write(22,"(i6, 3x, f13.8, 3x, f13.8)") aa, distPsi(aa, 1, ind1), distPsi(aa, 2, ind1)
	end do
end do

do ind1=1, rangodist(2)
	write(22,"(a6, 3x, f13.8, 3x, f13.8)") glob, distPsi(to+1, 1, ind1), distPsi(to+1, 2, ind1)
end do

!T
do aa=from, to
	do ind1=1, rangodist(3)
		write(23,"(i6, 3x, f13.8, 3x, f13.8)") aa, distT(aa, 1, ind1), distT(aa, 2, ind1)
	end do
end do

do ind1=1, rangodist(3)
	write(23,"(a6, 3x, f13.8, 3x, f13.8)") glob, distT(to+1, 1, ind1), distT(to+1, 2, ind1)
end do

!Ts
do aa=from, to
	do ind1=1, rangodist(4)
		write(24,"(i6, 3x, f13.8, 3x, f13.8)") aa, distTs(aa, 1, ind1), distTs(aa, 2, ind1)
	end do
end do

do ind1=1, rangodist(4)
	write(24,"(a6, 3x, f13.8, 3x, f13.8)") glob, distTs(to+1, 1, ind1), distTs(to+1, 2, ind1)
end do

!TCuad
do aa=from, to
	do ind1=1, rangodist(5)
		write(25,"(i6, 3x, f13.8, 3x, f13.8)") aa, distTCuad(aa, 1, ind1), distTCuad(aa, 2, ind1)
	end do
end do

do ind1=1, rangodist(5)
	write(25,"(a6, 3x, f13.8, 3x, f13.8)") glob, distTCuad(to+1, 1, ind1), distTCuad(to+1, 2, ind1)
end do

!Suma
do aa=from, to
	do ind1=1, rangodist(6)
		write(26,"(i6, 3x, f13.8, 3x, f13.8)") aa, distsuma(aa, 1, ind1), distsuma(aa, 2, ind1)
	end do
end do

do ind1=1, rangodist(6)
	write(26,"(a6, 3x, f13.8, 3x, f13.8)") glob, distsuma(to+1, 1, ind1), distsuma(to+1, 2, ind1)
end do

!Resta
do aa=from, to
	do ind1=1, rangodist(7)
		write(27,"(i6, 3x, f13.8, 3x, f13.8)") aa, distresta(aa, 1, ind1), distresta(aa, 2, ind1)
	end do
end do

do ind1=1, rangodist(7)
	write(27,"(a6, 3x, f13.8, 3x, f13.8)") glob, distresta(to+1, 1, ind1), distresta(to+1, 2, ind1)
end do

!TCuadsigno
do aa=from, to
	do ind1=1, rangodist(8)
		write(28,"(i6, 3x, f13.8, 3x, f13.8)") aa, distTCuadsigno(aa, 1, ind1), distTCuadsigno(aa, 2, ind1)
	end do
end do

do ind1=1, rangodist(8)
	write(28,"(a6, 3x, f13.8, 3x, f13.8)") glob, distTCuadsigno(to+1, 1, ind1), distTCuadsigno(to+1, 2, ind1)
end do


! Media y SD

do aa=from, to
	write(30, "(i6, 3x, a15, 3x, f15.8, 3x, f13.8)") aa, "dPhi", medias(aa,1), G(aa,1)
	write(30, "(i6, 3x, a15, 3x, f15.8, 3x, f13.8)") aa, "dPsi", medias(aa,2), G(aa,2)
	write(30, "(i6, 3x, a15, 3x, f15.8, 3x, f13.8)") aa, "dT", medias(aa,3), G(aa,3)
	write(30, "(i6, 3x, a15, 3x, f15.8, 3x, f13.8)") aa, "dT_por_signo", medias(aa,4), G(aa,4)
	write(30, "(i6, 3x, a15, 3x, f15.8, 3x, f13.8)") aa, "dT^2", medias(aa,5), G(aa,5)
	write(30, "(i6, 3x, a15, 3x, f15.8, 3x, f13.8)") aa, "dPhi+dPsi", medias(aa,6), G(aa,6)
	write(30, "(i6, 3x, a15, 3x, f15.8, 3x, f13.8)") aa, "dPhi-dPsi", medias(aa,7), G(aa,7)
	write(30, "(i6, 3x, a15, 3x, f15.8, 3x, f13.8)") aa, "dT_mismo_signo", mediaTsigno(aa,1), Gsigno(aa,1)
	write(30, "(i6, 3x, a15, 3x, f15.8, 3x, f13.8)") aa, "dT_dist_signo", mediaTsigno(aa,2), Gsigno(aa,2)
	write(30, "(i6, 3x, a15, 3x, f15.8, 3x, f13.8)") aa, "dT2_mismo_signo", mediaTCuadsigno(aa,1), GCuadsigno(aa,1)
	write(30, "(i6, 3x, a15, 3x, f15.8, 3x, f13.8)") aa, "dT2_dist_signo", mediaTCuadsigno(aa,2), GCuadsigno(aa,2)
end do

write(30, "(a6, 3x, a15, 3x, f15.8, 3x, f13.8)") glob, "dPhi", medias(to+1,1), G(to+1,1)
write(30, "(a6, 3x, a15, 3x, f15.8, 3x, f13.8)") glob, "dPsi", medias(to+1,2), G(to+1,2)
write(30, "(a6, 3x, a15, 3x, f15.8, 3x, f13.8)") glob, "dT", medias(to+1,3), G(to+1,3)
write(30, "(a6, 3x, a15, 3x, f15.8, 3x, f13.8)") glob, "dT por signo", medias(to+1,4), G(to+1,4)
write(30, "(a6, 3x, a15, 3x, f15.8, 3x, f13.8)") glob, "dT^2", medias(to+1,5), G(to+1,5)
write(30, "(a6, 3x, a15, 3x, f15.8, 3x, f13.8)") glob, "dPhi+dPsi", medias(to+1,6), G(to+1,6)
write(30, "(a6, 3x, a15, 3x, f15.8, 3x, f13.8)") glob, "dPhi-dPsi", medias(to+1,7), G(to+1,7)
write(30,"(a6,3x,a15,3x,f15.8,3x,f13.8)") glob,"dT_mismo_signo",mediaTsigno(to+1,1),Gsigno(to+1,1)
write(30,"(a6,3x,a15,3x,f15.8,3x,f13.8)") glob,"dT_dist_signo",mediaTsigno(to+1,2),Gsigno(to+1,2)
write(30,"(a6,3x,a15,3x,f15.8,3x,f13.8)") glob,"dT2_mismo_signo",mediaTCuadsigno(to+1,1),GCuadsigno(to+1,1)
write(30,"(a6,3x,a15,3x,f15.8,3x,f13.8)") glob,"dT2_dist_signo",mediaTCuadsigno(to+1,2),GCuadsigno(to+1,2)

end program distribuciones
