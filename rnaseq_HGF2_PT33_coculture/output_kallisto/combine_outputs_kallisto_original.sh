#!/bin/bash

# Directorio de salida para el archivo combinado
OUTPUT_DIR="/home/alexis/Desktop/github_script_finales_alexis/correr_kalisto/output_kallisto"

# Mover al directorio de salida
cd $OUTPUT_DIR

# Recorrer todas las carpetas en el directorio de salida
for FOLDER in */; do
    # Renombrar el archivo abundance.tsv al nombre de la carpeta
    mv "$FOLDER/abundance.tsv" "$FOLDER/${FOLDER%/}.tsv"
done

# Crear un archivo combinado con las columnas target_id y est_counts
echo -e "target_id\t"$(ls -1 */*.tsv | sed 's/\/.*//' | tr '\n' '\t') > combined_counts.tsv

# Combinar las columnas target_id y est_counts de todas las tablas abundance.tsv
cut -f 1 $(ls -1 */*.tsv | sed 's/\(.*\)/-f 2 \1/') | paste -d '\t' - >> combined_counts.tsv

