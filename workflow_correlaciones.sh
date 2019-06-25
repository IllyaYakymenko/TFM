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

direc_lista=""

#####################
while read linea ; do
#####################

#  Obtención de input para correlaciones.exe

cd auxdir
> ../correlaciones.in

##################################################
#
#  Fila 1: Número de ficheros
#  Fila 2: Lista de ficheros a utilizar
#  Fila 3: Número de filas de un fichero
#  Fila 4: Número de columnas de un fichero
#  Fila 5: PhiMin, PhiMax, PsiMin, PsiMax
#  Fila 6: Número de 
#  Fila 6: Phi, Psi de análisis
#  Fila 7: Resolución de mapa de correlaciones
#
##################################################

PhiRel="$(echo $linea | awk '{print $1}')"
PsiRel="$(echo $linea | awk '{print $2}')"

# Lista de ficheros (_raman de cada tray)
ficheros="$(ls *)"
echo $ficheros | wc -w >> ../correlaciones.in
echo $ficheros >> ../correlaciones.in

# Número de filas
fich1="$(echo $ficheros | awk '{print $1}')"
numfilas="$(echo $(wc -l $fich1))"
echo $numfilas| awk '{print $1}' >> ../correlaciones.in

# Número de columnas
head -n1 $fich1 | wc -w >> ../correlaciones.in

# PhiMin, PhiMax, PsiMin, PsiMax
# Límites del diagrama (-180, 180);(-180, 180)
limang="-180, 180, -180, 180"
echo $limang >> ../correlaciones.in

# Ángulos Phi y Psi de análisis respecto al aminoácido dado

echo "$PhiRel, $PsiRel" >> ../correlaciones.in

# Resolución del diagrama
echo $4 >> ../correlaciones.in

cp ../correlaciones.exe ./correlaciones.exe
cp ../correlaciones.in ./correlaciones.in
./correlaciones.exe < correlaciones.in

mv fort.20 "correlaciones_$1_$2_$3_resol$4_Phi$PhiRel-Psi$PsiRel-global.matrix"

numcol="$(head -n1 $fich1 | wc -w)"
numaa="$(expr $numcol / 2)"

from=1
to=$numaa

if [[ $PhiRel -gt 0 || $PsiRel -gt 0 ]] ; then
	if [[ $PhiRel -gt $PsiRel ]] ; then
		to=$(($numaa-$PhiRel))
	else
		to=$(($numaa-$PsiRel))
	fi
fi

if [[ $PhiRel -lt 0 || $PsiRel -lt 0 ]] ; then
	if [[ $PhiRel -lt $PsiRel ]] ; then
		from=$((1-$PhiRel))
	else
		from=$((1-$PsiRel))
	fi
fi

for aa in $(seq $from $to) ; do	
	ch=$(($aa+20))
	mv "fort.$ch" "correlaciones_$1_$2_$3_resol$4_Phi$PhiRel-Psi$PsiRel-aa_$aa.matrix"
done

rm -rf ../correlaciones_$1_$2_$3_resol$4_Phi$PhiRel-Psi$PsiRel.matrices
mkdir ../correlaciones_$1_$2_$3_resol$4_Phi$PhiRel-Psi$PsiRel.matrices
cp *.matrix ../correlaciones_$1_$2_$3_resol$4_Phi$PhiRel-Psi$PsiRel.matrices

direc="correlaciones_$1_$2_$3_resol$4_Phi$PhiRel-Psi$PsiRel.matrices"

direc_lista="$direc_lista $direc"

rm correlaciones.exe
rm correlaciones.in
rm *.matrix

cd ..

#############################
done < angulosrelativos.input
#############################

mkdir correlaciones_$1_$2_$3_resol$4 2>/dev/null

for directorio in $direc_lista ; do
	mv $directorio ./correlaciones_$1_$2_$3_resol$4/$directorio
done

#rm -rf auxdir













