#!/bin/bash

################################################################
#  
#  $1 = Nombre archivos
#  $2 = Numero inicio ficheros
#  $3 = Numero fin ficheros
#  $4 = Resolución
#  $5 = PhiRel
#  $6 = PsiRel
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

#  Obtención de input para mapa_modulos_totales.exe

cd auxdir
> ../mapa_modulos_totales.in

##############################################
#
#  Fila 1: Número de ficheros
#  Fila 2: Lista de ficheros a utilizar
#  Fila 3: Número de filas de un fichero
#  Fila 4: Número de columnas de un fichero
#  Fila 5: PhiMin, PhiMax, PsiMin, PsiMax
#  Fila 6: Resolución del mapa
#  Fila 7: PhiRel, PsiRel
#
##############################################


# Lista de ficheros (_raman de cada tray)
ficheros="$(ls *)"
echo $ficheros | wc -w >> ../mapa_modulos_totales.in
echo $ficheros >> ../mapa_modulos_totales.in

# Número de filas
fich1="$(echo $ficheros | awk '{print $1}')"
numfilas="$(echo $(wc -l $fich1))"
echo $numfilas| awk '{print $1}' >> ../mapa_modulos_totales.in

# Número de columnas
head -n1 $fich1 | wc -w >> ../mapa_modulos_totales.in

# PhiMin, PhiMax, PsiMin, PsiMax
# Límites del diagrama (-180, 180);(-180, 180)
limang="-180, 180, -180, 180"
echo $limang >> ../mapa_modulos_totales.in

# Resolución del mapa
echo $4 >> ../mapa_modulos_totales.in

#PhiRel, PsiRel
echo $5, $6 >> ../mapa_modulos_totales.in

cp ../mapa_modulos_totales.exe ./mapa_modulos_totales.exe
cp ../mapa_modulos_totales.in ./mapa_modulos_totales.in
./mapa_modulos_totales.exe < mapa_modulos_totales.in

mv fort.20 "mapa_totales_$1_$2_$3_resol$4_$5_$6_global.matrix"

numcol="$(head -n1 $fich1 | wc -w)"
numaa="$(expr $numcol / 2)" 

for aa in $(seq 1 $numaa) ; do	
	ch=$(($aa+20))
	mv "fort.$ch" "mapa_totales_$1_$2_$3_resol$4_$5_$6_aa_$aa.matrix" 2> /dev/null
done

mv fort.50 "mapa_totalesCuadrado_$1_$2_$3_resol$4_$5_$6_global.matrix"

for aa in $(seq 1 $numaa) ; do	
	ch=$(($aa+50))
	mv "fort.$ch" "mapa_totalesCuadrado_$1_$2_$3_resol$4_$5_$6_aa_$aa.matrix" 2> /dev/null
done

# Las salidas 20->30 generan los mapas de módulos D y las salidas 50->60 los de módulos D^2

rm -rf ../mapa_totales_$1_$2_$3_resol$4_$5_$6_matrices
mkdir ../mapa_totales_$1_$2_$3_resol$4_$5_$6_matrices
cp *.matrix ../mapa_totales_$1_$2_$3_resol$4_$5_$6_matrices

cd ..
rm -rf auxdir













