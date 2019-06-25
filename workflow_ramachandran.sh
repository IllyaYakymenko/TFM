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

#  Obtención de input para ramachandran.exe

cd auxdir
> ../ramachandran.in

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
echo $ficheros | wc -w >> ../ramachandran.in
echo $ficheros >> ../ramachandran.in

# Número de filas
fich1="$(echo $ficheros | awk '{print $1}')"
numfilas="$(echo $(wc -l $fich1))"
echo $numfilas| awk '{print $1}' >> ../ramachandran.in

# Número de columnas
head -n1 $fich1 | wc -w >> ../ramachandran.in

# PhiMin, PhiMax, PsiMin, PsiMax
# Límites del diagrama (-180, 180);(-180, 180)
limang="-180, 180, -180, 180"
echo $limang >> ../ramachandran.in

# Resolución del diagrama
echo $4 >> ../ramachandran.in

cp ../ramachandran.exe ./ramachandran.exe
cp ../ramachandran.in ./ramachandran.in
./ramachandran.exe < ramachandran.in

mv fort.20 "ramachandran_$1_$2_$3_resol$4_global.matrix"

numcol="$(head -n1 $fich1 | wc -w)"
numaa="$(expr $numcol / 2)" 

for aa in $(seq 1 $numaa) ; do	
	ch=$(($aa+20))
	mv "fort.$ch" "ramachandran_$1_$2_$3_resol$4_aa_$aa.matrix"
done

rm -rf ../ramachandran_$1_$2_$3_resol$4_matrices
mkdir ../ramachandran_$1_$2_$3_resol$4_matrices
cp *.matrix ../ramachandran_$1_$2_$3_resol$4_matrices

cd ..
rm -rf auxdir













