#!/bin/bash

################################################################
#  
#  $1 = Nombre archivos
#  $2 = Numero inicio ficheros
#  $3 = Numero fin ficheros
#  $4 = Resolución
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

#  Obtención de input para mapasDE.exe

cd auxdir
> ../mapasDE.in

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
echo $ficheros | wc -w >> ../mapasDE.in
echo $ficheros >> ../mapasDE.in

# Número de filas
fich1="$(echo $ficheros | awk '{print $1}')"
numfilas="$(echo $(wc -l $fich1))"
echo $numfilas| awk '{print $1}' >> ../mapasDE.in

# Número de columnas
head -n1 $fich1 | wc -w >> ../mapasDE.in

# Resolución del diagrama
echo $4 >> ../mapasDE.in

cp ../mapasDE.exe ./mapasDE.exe
cp ../mapasDE.in ./mapasDE.in
./mapasDE.exe < mapasDE.in

mv fort.20 "Mapa_Sigma_dPhi(i)_$1_$2_$3_resol$4.matrix"
mv fort.21 "Mapa_Sigma_dPsi(i)_$1_$2_$3_resol$4.matrix"
mv fort.22 "Mapa_Sigma_dPhi(i+1)_$1_$2_$3_resol$4.matrix"
mv fort.23 "Mapa_Sigma_dPsi(i-1)_$1_$2_$3_resol$4.matrix"

rm -rf ../mapasDE_$1_$2_$3_resol$4_matrices
mkdir ../mapasDE_$1_$2_$3_resol$4_matrices
cp *.matrix ../mapasDE_$1_$2_$3_resol$4_matrices

cd ..
rm -rf auxdir













