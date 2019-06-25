#!/bin/bash

################################################################
#  
#  $1 = Nombre archivos
#  $2 = Numero inicio ficheros
#  $3 = Numero fin ficheros
#  $4 = Intervalo de registro
#  $5 = Rango límite
#  $6 = Phi centro
#  $7 = Psi centro
#  $8 = Rango Phi
#  $9 = Rango Psi
#  $10 = RelPhi
#  $11 = RelPsi
#
################################################################

rm -rf auxdir
mkdir auxdir

#  Copia de los archivos raman.dat

padding=5
zi=""
for ((a=1; a <= "$padding" ; a++))
do
zi=$zi"0"
done
lista=""
for ((a="$2"; a <= "$3" ; a++))
do
la=${#a}
let "cla=$padding-$la"
z=$zi
b=${z:$la:$cla}$a
lista=$lista$b" "
done

actualdir="$(pwd)"

maximo="$(($3-$2+1))"
cont_fich=0
echo "Descomprimiendo [0/$maximo]"

for a in $lista; do
    cont_fich="$(($cont_fich+1))"
    arch=$1e$a\_raman.dat.gz
    dirtray="$actualdir/tray0/tray$a"
    cp $dirtray/$arch ./auxdir
    gunzip ./auxdir/$arch
    echo "Descomprimiendo [$cont_fich/$maximo]"
done

#  Obtención de input para ramachandran.exe

cd auxdir
> ../dens_desplazamientos.in

##############################################
#
#  Fila 1: Número de ficheros
#  Fila 2: Lista de ficheros a utilizar
#  Fila 3: Número de filas de un fichero
#  Fila 4: Número de columnas de un fichero
#  Fila 5: PhiMin, PhiMax, PsiMin, PsiMax
#  Fila 6: Resolución Ramachandran
#
##############################################


# Lista de ficheros (_raman de cada tray)
ficheros="$(ls *)"
echo $ficheros | wc -w >> ../dens_desplazamientos.in
echo $ficheros >> ../dens_desplazamientos.in

# Número de filas
fich1="$(echo $ficheros | awk '{print $1}')"
numfilas="$(echo $(wc -l $fich1))"
echo $numfilas| awk '{print $1}' >> ../dens_desplazamientos.in

# Número de columnas
head -n1 $fich1 | wc -w >> ../dens_desplazamientos.in

# PhiMin, PhiMax, PsiMin, PsiMax
# Límites del diagrama (-180, 180);(-180, 180)
limang="$(echo $6, $7, $8, $9)"
echo $limang >> ../dens_desplazamientos.in

# Intervalo de registro
echo $4 >> ../dens_desplazamientos.in

# Rango límite
echo $5 >> ../dens_desplazamientos.in

# PhiRel, PsiRel
echo ${10}, ${11} >> ../dens_desplazamientos.in

cp ../dens_desplazamientos.exe ./dens_desplazamientos.exe
cp ../dens_desplazamientos.in ./dens_desplazamientos.in
./dens_desplazamientos.exe < dens_desplazamientos.in

mv fort.20 "densPhiPsi_$1_$2_$3_resol$4_lim$5_Phi$6_$8_${10}_Psi$7_$9_${11}.matrix"
mv fort.30 "denssumaresta_$1_$2_$3_resol$4_lim$5_Phi$6_$8_${10}_Psi$7_$9_${11}.matrix"

rm -rf ../dens_desplazamientos_$1_$2_$3_resol$4_lim$5_Phi$6_$8_${10}_Psi$7_$9_${11}
mkdir ../dens_desplazamientos_$1_$2_$3_resol$4_lim$5_Phi$6_$8_${10}_Psi$7_$9_${11}
cp *.matrix ../dens_desplazamientos_$1_$2_$3_resol$4_lim$5_Phi$6_$8_${10}_Psi$7_$9_${11}

cd ..
rm -rf auxdir













