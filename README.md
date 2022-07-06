# Monkeypox_assembly
collection of codes for monkeypox genome assembly

## STEP 0: just in case generate a naked script.sh file
```r
cd $HOME && touch script.sh && chmod +x script.sh ;
cd $HOME && echo '#!/bin/bash' > script.sh && echo '# -*- ENCODING: UTF-8 -*-' >> script.sh ;
mv script.sh paste/your/working/directory/ ;
```

## STEP 1: copy and paste the following code in the .sh file adding "exit" in the last line
```r
#1#indexado y alineado con el genoma de referencia#
bwa index NC_063383.1.fasta ; 
#ojo: para genomas codificados como CODIGO_INS_1.fastq.gz y CODIGO_INS_2.fastq.gz#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _1.fastq.gz)
r2=${prefix}_2.fastq.gz

#2# instrucciones para generar el archivo .bam#
bwa mem -t 15 NC_063383.1.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 15 -bS -T NC_063383.1.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 15 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 15 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 15 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 15 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 15 ${prefix}.bam ;

#3#remocion de intermedios#
rm *_uno.bam *_uno.sam *_unoa.bam *_dosa.bam *_tresa.bam *_cuatroa.bam ;

#4#obtencion del genoma consenso#
samtools mpileup -aa -A -d 0 -Q 0 ${prefix}.bam | ivar consensus -p ${prefix}.fasta -q 25 -t 0.6 -m 20 ; 

#5#obtencion de archivos bam con solo mapeados#
samtools view -b -F 4 ${prefix}.bam > ${prefix}.mapped.bam ; 
samtools fastq -1 ${prefix}_f.fq -2 ${prefix}_r.fq -0 /dev/null -s /dev/null -n ${prefix}.mapped.bam ; 
rm *.fastq.gz.fa *.fastq.gz.mapped.bam *.fastq.gz.qual.txt *.fastq.gz_f.fq *.fastq.gz_r.fq ;
samtools index -@ 15 ${prefix}.mapped.bam ; 
done ;

#6#mover los archivos#
mkdir mapped ; 
mv *.mapped.bam *.fa *.fq *.mapped.bam.bai mapped/ ;
cat mapped/*.fa > mapped/genomes.fasta ; 
aliview mapped/genomes.fasta ;
```

## STEP 2: count Ns for every genome
```r
#1#set names#
for r1 in *fa
do
prefix=$(basename $r1 .fa)
#2#estimate Ns#
seqtk comp $r1 | awk '{x=$3+$4+$5+$6;y=$2;print $1,y-x,y,(y-x)/y}' > ${prefix}_count.txt ;
done ; 
#3#merge data#
cat *_count.txt > n_count.txt ;
```

## STEP 3: qualimap report
download qualimap from http://qualimap.conesalab.org/, uncompress and paste in working directory
```r
#1#set names#
for r1 in *mapped.bam
do
prefix=$(basename $r1 .mapped.bam)

#2# correr qualimap para cada *.bam#
./qualimap bamqc -bam $r1 -outfile ${prefix}.pdf ; 
mv ${prefix}.mapped_stats/${prefix}.pdf . ;
done ; 
```

## STEP 4: repeat ivar with m = 5
```r
#1#set names#
for r1 in *bam
do
prefix=$(basename $r1 .bam)
#2#estimate Ns#
samtools mpileup -aa -A -d 0 -Q 0 $r1 | ivar consensus -p ${prefix}.fasta -q 25 -t 0.6 -m 5 ;
done ; 
## Repeat STEP 2
```
