#!/bin/bash

################################################################
#  
#  $1 = Nombre archivos
#  $2 = Numero inicio ficheros
#  $3 = Numero fin ficheros
#  $4 = Registro de datos desde $4% del valor máximo/mínimo
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

cd auxdir

cp ../mayores_modulos.exe ./mayores_modulos.exe
> mayores_modulos.in

######################################################
#
#  Fila 1: Número de ficheros
#  Fila 2: Lista de ficheros
#  Fila 3: Número de filas
#  Fila 4: Número de columnas
#  Fila 5: % de datos extremos
#
######################################################


# Lista de ficheros (_raman de cada tray)
ficheros="$(ls *.dat)"
echo $ficheros | wc -w >> mayores_modulos.in
echo $ficheros >> mayores_modulos.in

# Número de filas
fich1="$(echo $ficheros | awk '{print $1}')"
numfilas="$(echo $(wc -l $fich1))"
echo $numfilas| awk '{print $1}' >> mayores_modulos.in

# Número de columnas
head -n1 $fich1 | wc -w >> mayores_modulos.in

# Significancia
echo $4 >> mayores_modulos.in

./mayores_modulos.exe < mayores_modulos.in

mv fort.20 ../"Mayores_modulos_$2_$3_registrado$4%.dat"
cd ..

rm -rf auxdir











