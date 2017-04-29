MIBLAB=/opt/miblab
OPENCV=/opt/opencv
LAPACK=/home/jphan/Projects/CodeTree/lapack-3.2.1
MYSQLINC=/usr/include/mysql
MYSQLLIB=/usr/lib64/mysql
BOOST=/home/jphan/Projects/CodeTree/boost_1_42_0

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
./obj/seqc_db.o: ./src/seqc_db.cpp ./include/seqc_db.h
	$(CC) -c -O2 -o ./obj/seqc_db.o \
		./src/seqc_db.cpp \
		-I$(MIBLAB)/include

#--------------
seqc_upload_quantification: ./obj/seqc_upload_quantification.o ./obj/seqc_db.o
	$(CC) -o ./seqc_upload_quantification ./obj/seqc_upload_quantification.o \
		./obj/seqc_db.o \
		-L$(MIBLAB)/lib -lmiblab \
		-L$(MYSQLLIB) -lmysqlclient

./obj/seqc_upload_quantification.o: ./src/seqc_upload_quantification.cpp
	$(CC) -c -O2 -o ./obj/seqc_upload_quantification.o \
		./src/seqc_upload_quantification.cpp \
		-I$(MIBLAB)/include \
		-I$(MYSQLINC)

#--------------
seqc_extract_expression: ./obj/seqc_extract_expression.o ./obj/seqc_db.o
	$(CC) -o ./seqc_extract_expression ./obj/seqc_extract_expression.o \
		./obj/seqc_db.o \
		-L$(MIBLAB)/lib -lmiblab \
		-L$(MYSQLLIB) -lmysqlclient

./obj/seqc_extract_expression.o: ./src/seqc_extract_expression.cpp
	$(CC) -c -O2 -o ./obj/seqc_extract_expression.o \
		./src/seqc_extract_expression.cpp \
		-I$(MIBLAB)/include \
		-I$(MYSQLINC)

#--------------
clean:
	rm $(APPS) -f
	rm ./obj/*.o -f
