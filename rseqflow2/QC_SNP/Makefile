all: Merge_genome_transcriptome_singleEnd Merge_genome_transcriptome_pairEnd


Merge_genome_transcriptome_singleEnd: Merge_genome_transcriptome_singleEnd.o 
	g++ -o Merge_genome_transcriptome_singleEnd Merge_genome_transcriptome_singleEnd.o -Wall -Wno-deprecated -g 

Merge_genome_transcriptome_singleEnd.o: Merge_genome_transcriptome_singleEnd.cpp Merge_genome_transcriptome.h
	g++ -c Merge_genome_transcriptome_singleEnd.cpp -Wall -Wno-deprecated -g  

Merge_genome_transcriptome_pairEnd: Merge_genome_transcriptome_pairEnd.o
	g++ -o Merge_genome_transcriptome_pairEnd Merge_genome_transcriptome_pairEnd.o -Wall -Wno-deprecated -g

Merge_genome_transcriptome_pairEnd.o: Merge_genome_transcriptome_pairEnd.cpp Merge_genome_transcriptome.h
	g++ -c Merge_genome_transcriptome_pairEnd.cpp -Wall -Wno-deprecated -g


clean: 
	rm -f Merge_genome_transcriptome_singleEnd Merge_genome_transcriptome_pairEnd

