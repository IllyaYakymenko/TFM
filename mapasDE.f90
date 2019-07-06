! #################################################################################
! #
! #   Obtención de mapas de desplazamientos estándar
! #   
! #   Input: Channel 10 (fichero _raman.dat)
! #   Output: Channel 20: Mapa Sigma Phi(i) 
! #	                  21: Mapa Sigma Psi(i)
! #                   22: Mapa Sigma Phi(i+1)
! #                   23: Mapa Sigma Psi(i-1)
! #
! #################################################################################

program mapasDE

implicit none

! Datos de entrada (archivo.in)

integer :: numfiles, resol, nrow, ncol
character*200, allocatable :: filenames(:)

! Definición de la matriz del diagrama

double precision, allocatable :: SumaPhi(:,:), SumaPsi(:,:), SumaPhiPost(:,:), SumaPsiAnt(:,:)
double precision, allocatable :: NPhiPsi(:,:), NPost(:,:), NAnt(:,:)
double precision, allocatable :: MediaPhi(:,:), MediaPsi(:,:), MediaPhiPost(:,:), MediaPsiAnt(:,:)
double precision, allocatable :: SigmaPhi(:,:), SigmaPsi(:,:), SigmaPsiAnt(:,:), SigmaPhiPost(:,:)
double precision, allocatable :: SSigmaPhi(:,:), SSigmaPsi(:,:), SSigmaPsiAnt(:,:), SSigmaPhiPost(:,:)
double precision :: rangoPhi, rangoPsi, dPhi, dPsi, dPhiPost, dPsiAnt
integer :: diagPhi, diagPsi, numaa
double precision, allocatable :: diedros(:), diedrosant(:), desps(:)

! Variables de control

integer ch, row, col, fileind, rowiter, aa, i, diagaa, di

read(5,*) numfiles

allocate(filenames(numfiles))

read(5,*) filenames
read(5,*) nrow
read(5,*) ncol
read(5,*) resol

numaa = ncol/2

allocate(diedros(ncol))
allocate(diedrosant(ncol))
allocate(desps(ncol))
allocate(SumaPhi(resol, resol))
allocate(SumaPsi(resol, resol))
allocate(SumaPhiPost(resol, resol))
allocate(SumaPsiAnt(resol, resol))
allocate(NPhiPsi(resol, resol))
allocate(NAnt(resol, resol))
allocate(NPost(resol, resol))
allocate(MediaPhi(resol, resol))
allocate(MediaPsi(resol, resol))
allocate(MediaPhiPost(resol, resol))
allocate(MediaPsiAnt(resol, resol))
allocate(SigmaPhi(resol, resol))
allocate(SigmaPsi(resol, resol))
allocate(SigmaPhiPost(resol, resol))
allocate(SigmaPsiAnt(resol, resol))
allocate(SSigmaPhi(resol, resol))
allocate(SSigmaPsi(resol, resol))
allocate(SSigmaPhiPost(resol, resol))
allocate(SSigmaPsiAnt(resol, resol))

! Matrices vacias 

do row=1, resol
	do col=1, resol
		SumaPhi(row,col) = 0.d0
		SumaPsi(row,col) = 0.d0
		SumaPhiPost(row,col) = 0.d0
		SumaPsiAnt(row,col) = 0.d0
		SSigmaPhi(row,col) = 0.d0
		SSigmaPsi(row,col) = 0.d0
		SSigmaPhiPost(row,col) = 0.d0
		SSigmaPsiAnt(row,col) = 0.d0
		NPhiPsi(row, col) = 0.d0
		NAnt(row, col) = 0.d0
		NPost(row, col) = 0.d0
	end do
end do

rangoPhi = 360/resol
rangoPsi = 360/resol

print*, 'Calculando medias'

! Registro de sumatorios para el cálculo de medias (primer bucle)

do fileind=1, numfiles
	open(10, file = filenames(fileind))
	read(10,*) diedrosant

	do rowiter=2, nrow
		read(10,*) diedros

		do di=1, ncol
		desps(di) = (diedros(di) - diedrosant(di))-1*int((diedros(di) - diedrosant(di))/180)*360	
		end do

		do aa=1, numaa
			! Posición en las matrices
			diagPhi = int((diedros(2*aa-1)+180)/rangoPhi)+1
			diagPsi = int((diedros(2*aa)+180)/rangoPsi)+1
			
			if(diagPhi.ge.resol) then
				diagPhi=resol
			end if

			if(diagPsi.ge.resol) then
				diagPsi=resol
			end if

			! Sumatorios de desplazamientos de ángulos pertenecientes al mismo residuo
			dPhi = desps(2*aa -1)
			dPsi = desps(2*aa)
			NPhiPsi(diagPsi, diagPhi) = NPhiPsi(diagPsi, diagPhi) +1
			SumaPhi(diagPsi, diagPhi) = SumaPhi(diagPsi, diagPhi) + dPhi
			SumaPsi(diagPsi, diagPhi) = SumaPsi(diagPsi, diagPhi) + dPsi

			! Sumatorios de desplazamientos de Psi del residuo anterior
			! El residuo 1 no posee ángulo Psi anterior al residuo 1
			if(aa.ge.2) then
				NAnt(diagPsi, diagPhi) = NAnt(diagPsi, diagPhi) +1
				dPsiAnt = desps(2*aa-2)
				SumaPsiAnt(diagPsi, diagPhi) = SumaPsiAnt(diagPsi, diagPhi) + dPsiAnt
			end if

			! Sumatorios de desplazamientos de Phi del residuo posterior
			! El residuo numaa no posee ángulo Phi anterior al residuo numaa
			if(aa.le.numaa-1) then
				NPost(diagPsi, diagPhi) = NPost(diagPsi, diagPhi) +1
				dPhiPost = desps(2*aa+1)
				SumaPhiPost(diagPsi, diagPhi) = SumaPhiPost(diagPsi, diagPhi) + dPhiPost
			end if
			
		end do
		diedrosant = diedros
	end do
	close(10)
	print '("Calculando medias: [", i3, "/", i3, "]")', fileind, numfiles
end do

! Cálculo de las medias de los desplazamientos

do row=1, resol
	do col=1, resol
		if(NPhiPsi(row, col).gt.0) then
			MediaPhi(row, col) = SumaPhi(row, col) / NPhiPsi(row, col)
			MediaPsi(row, col) = SumaPsi(row, col) / NPhiPsi(row, col)
		else
			MediaPhi(row, col) = 0.d0
			MediaPsi(row, col) = 0.d0
		end if

		if(NPost(row, col).gt.0) then
			MediaPhiPost(row, col) = SumaPhiPost(row, col) / NPost(row, col)
		else
			MediaPhiPost(row, col) = 0.d0
		end if

		if(NAnt(row, col).gt.0) then
			MediaPsiAnt(row, col) = SumaPsiAnt(row, col) / NAnt(row, col)
		else
			MediaPsiAnt(row, col) = 0.d0
		end if

	end do
end do

print*, 'Registrando desviación estándar'

! Registro de sumatorios para el cálculo de desviaciones estándar (segundo bucle)

do fileind=1, numfiles
	open(10, file = filenames(fileind))
	read(10,*) diedrosant

	do rowiter=2, nrow
		read(10,*) diedros

		do di=1, ncol
		desps(di) = (diedros(di) - diedrosant(di))-1*int((diedros(di) - diedrosant(di))/180)*360	
		end do

		do aa=1, numaa
			! Posición en las matrices
			diagPhi = int((diedros(2*aa-1)+180)/rangoPhi)+1
			diagPsi = int((diedros(2*aa)+180)/rangoPsi)+1
			
			if(diagPhi.ge.resol) then
				diagPhi=resol
			end if

			if(diagPsi.ge.resol) then
				diagPsi=resol
			end if

			! Ángulos del mismo residuo
			dPhi = desps(2*aa -1)
			dPsi = desps(2*aa)
			SSigmaPhi(diagPsi,diagPhi)=SSigmaPhi(diagPsi,diagPhi)+(MediaPhi(diagPsi,diagPhi)-dPhi)**2
			SSigmaPsi(diagPsi,diagPhi)=SSigmaPsi(diagPsi,diagPhi)+(MediaPsi(diagPsi,diagPhi)-dPsi)**2

			! Ángulo del residuo anterior
			if(aa.ge.2) then
				dPsiAnt = desps(2*aa-2)
				SSigmaPsiAnt(diagPsi,diagPhi)=SSigmaPsiAnt(diagPsi,diagPhi)+(MediaPsiAnt(diagPsi,diagPhi)-dPsiAnt)**2
			end if

			! Ángulo del residuo posterior
			if(aa.le.numaa-1) then
				dPhiPost = desps(2*aa+1)
				SSigmaPhiPost(diagPsi,diagPhi)=SSigmaPhiPost(diagPsi,diagPhi)+(MediaPhiPost(diagPsi,diagPhi)-dPhiPost)**2
			end if
			
		end do
		diedrosant = diedros
	end do
	close(10)
	print '("Registrando desviación estándar: [", i3, "/", i3, "]")', fileind, numfiles
end do

! Calculando las desviaciones estándar

do row=1, resol
	do col=1, resol
		if(NPhiPsi(row, col).gt.0) then
			SigmaPhi(row, col) = sqrt(SSigmaPhi(row, col)/NPhiPsi(row, col))
			SigmaPsi(row, col) = sqrt(SSigmaPsi(row, col)/NPhiPsi(row, col))
		else
			SigmaPhi(row, col) = -10
			SigmaPsi(row, col) = -10
		end if

		if(NPost(row, col).gt.0) then
			SigmaPhiPost(row, col) = sqrt(SSigmaPhiPost(row, col)/NPost(row, col))
		else
			SigmaPhiPost(row, col) = -10
		end if

		if(NAnt(row, col).gt.0) then
			SigmaPsiAnt(row, col) = sqrt(SSigmaPsiAnt(row, col)/NAnt(row, col))
		else
			SigmaPsiAnt(row, col) = -10
		end if
	end do
end do

! Escritura los resultados

do row=1, resol
	write(20,*) (SigmaPhi(row, col), col=1, resol)
	write(21,*) (SigmaPsi(row, col), col=1, resol)
	write(22,*) (SigmaPhiPost(row, col), col=1, resol)
	write(23,*) (SigmaPsiAnt(row, col), col=1, resol)
end do

end program mapasDE




