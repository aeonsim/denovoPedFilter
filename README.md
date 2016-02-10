denovoPedFilter
===============

A multi-generation based pedigree filter for the identification and selection of de novo mutations.

Uses mendelian phasing in 3 & 4 generation pedigrees & outputs various metrics about possible de novos. 

While designed for running on population VCF files containing larger pedigrees (ie parents, proband and children or grandparents, parents, proband & children), it will also handle the simple case of a basic Trios though the resulting list of denovo mutations will have a much higher false positive rate. Currently supports VCFs from Freebayes, GATK & Platypus.

Download the JAR from here: https://github.com/aeonsim/denovoPedFilter/blob/master/target/scala-2.11/pedigreeFilter.jar
