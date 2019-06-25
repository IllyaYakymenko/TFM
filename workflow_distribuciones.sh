#!/bin/bash

################################################################
#  
#  $1 = Nombre archivos
#  $2 = Numero inicio ficheros
#  $3 = Numero fin ficheros
#  $4 = Intervalo de registro
#  $5 = Phi centro
#  $6 = Psi centro
#  $7 = Rango Phi
#  $8 = Rango Psi
#  $9 = RelPhi
#  $10 = RelPsi
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
> ../distribuciones.in

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
echo $ficheros | wc -w >> ../distribuciones.in
echo $ficheros >> ../distribuciones.in

# Número de filas
fich1="$(echo $ficheros | awk '{print $1}')"
numfilas="$(echo $(wc -l $fich1))"
echo $numfilas| awk '{print $1}' >> ../distribuciones.in

# Número de columnas
head -n1 $fich1 | wc -w >> ../distribuciones.in

# PhiMin, PhiMax, PsiMin, PsiMax
# Límites del diagrama (-180, 180);(-180, 180)
limang="$(echo $5, $6, $7, $8)"
echo $limang >> ../distribuciones.in

# Intervalo de registro
echo $4 >> ../distribuciones.in

# RelPhi, RelPsi
echo $9, ${10} >> ../distribuciones.in

cp ../distribuciones.exe ./distribuciones.exe
cp ../distribuciones.in ./distribuciones.in
./distribuciones.exe < distribuciones.in

mv fort.21 "distdPhi_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}.hist"
mv fort.22 "distdPsi_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}.hist"
mv fort.23 "distdT_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}.hist"
mv fort.24 "distdTsigno_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}.hist"
mv fort.25 "distdTCuadrado_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}.hist"
mv fort.26 "distdPhi+dPsi_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}.hist"
mv fort.27 "distdPhi-dPsi_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}.hist"
mv fort.28 "distdT2signo_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}.hist"

mv fort.30 "medias_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}.dat"

rm -rf ../distribuciones_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}
mkdir ../distribuciones_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}
cp *.hist ../distribuciones_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}
cp medias*.dat ../distribuciones_$1_$2_$3_resol$4_Phi$5_$7_$9_Psi$6_$8_${10}

cd ..
rm -rf auxdir













