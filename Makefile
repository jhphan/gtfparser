MIBLAB=/opt/miblab

CC=g++

APPS=	seqc_upload_quantification \
	seqc_extract_expression \
	transcriptome_to_genome \
	transcriptome_to_genome2 \
	read_sam_flags \
	filter_sam2 \
	merge_exons \
	compute_exon_gene_length \
	filter_duplicates

all: $(APPS)


./obj/sam.o: ./src/sam.cpp
	$(CC) -c -O2 -o ./obj/sam.o \
		./src/sam.cpp \
		-I$(MIBLAB)/include

./obj/transcript.o: ./src/transcript.cpp
	$(CC) -c -O2 -o ./obj/transcript.o \
		./src/transcript.cpp \
		-I$(MIBLAB)/include

#--------------

./obj/filter_sam2.o: ./src/filter_sam2.cpp
	$(CC) -c -O2 -o ./obj/filter_sam2.o \
		./src/filter_sam2.cpp \
		-I$(MIBLAB)/include

filter_sam2: ./obj/filter_sam2.o
	$(CC) -o ./filter_sam2 \
		./obj/filter_sam2.o \
		-L$(MIBLAB)/lib -lmiblab

#--------------

./obj/filter_duplicates.o: ./src/filter_duplicates.cpp
	$(CC) -c -O2 -o ./obj/filter_duplicates.o \
		./src/filter_duplicates.cpp \
		-I$(MIBLAB)/include

filter_duplicates: ./obj/filter_duplicates.o
	$(CC) -o ./filter_duplicates \
		./obj/filter_duplicates.o \
		-L$(MIBLAB)/lib -lmiblab

# --------------

./obj/read_sam_flags.o: ./src/read_sam_flags.cpp
	$(CC) -c -O2 -o ./obj/read_sam_flags.o \
		./src/read_sam_flags.cpp \
		-I$(MIBLAB)/include

read_sam_flags: ./obj/read_sam_flags.o
	$(CC) -o ./read_sam_flags \
		./obj/read_sam_flags.o \
		-L$(MIBLAB)/lib -lmiblab

# --------------
compute_exon_gene_length: ./obj/compute_exon_gene_length.o ./obj/transcript.o
	$(CC) -o ./compute_exon_gene_length \
		./obj/compute_exon_gene_length.o \
		./obj/transcript.o \
		-L$(MIBLAB)/lib -lmiblab

./obj/compute_exon_gene_length.o: ./src/compute_exon_gene_length.cpp
	$(CC) -c -O2 -o ./obj/compute_exon_gene_length.o \
		./src/compute_exon_gene_length.cpp \
		-I$(MIBLAB)/include 

# --------------
merge_exons: ./obj/merge_exons.o ./obj/transcript.o
	$(CC) -o ./merge_exons \
		./obj/merge_exons.o \
		./obj/transcript.o \
		-L$(MIBLAB)/lib -lmiblab

./obj/merge_exons.o: ./src/merge_exons.cpp
	$(CC) -c -O2 -o ./obj/merge_exons.o \
		./src/merge_exons.cpp \
		-I$(MIBLAB)/include 

# --------------
transcriptome_to_genome: ./obj/transcriptome_to_genome.o ./obj/transcript.o ./obj/sam.o
	$(CC) -o ./transcriptome_to_genome \
		./obj/transcriptome_to_genome.o \
		./obj/sam.o \
		./obj/transcript.o \
		-L$(MIBLAB)/lib -lmiblab

./obj/transcriptome_to_genome.o: ./src/transcriptome_to_genome.cpp
	$(CC) -c -O2 -o ./obj/transcriptome_to_genome.o \
		./src/transcriptome_to_genome.cpp \
		-I$(MIBLAB)/include 

# --------------
transcriptome_to_genome2: ./obj/transcriptome_to_genome2.o ./obj/transcript.o ./obj/sam.o
	$(CC) -o ./transcriptome_to_genome2 \
		./obj/transcriptome_to_genome2.o \
		./obj/sam.o \
		./obj/transcript.o \
		-L$(MIBLAB)/lib -lmiblab

./obj/transcriptome_to_genome2.o: ./src/transcriptome_to_genome2.cpp
	$(CC) -c -O2 -o ./obj/transcriptome_to_genome2.o \
		./src/transcriptome_to_genome2.cpp \
		-I$(MIBLAB)/include

#--------------
clean:
	rm $(APPS) -f
	rm ./obj/*.o -f
