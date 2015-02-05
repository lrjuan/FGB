##############################

Local Installation Instruction of Family Genome Browser

##############################

The Family Genome Browser (FGB) can be installed in Linux/Mac OSX environment. 

1. JAVA 6 and Tomcat 6 or later versions should be ready on the server.

Here are installation instruction for

JAVA:
http://www.java.com/en/download/help/index_installing.xml

Tomcat:
https://tomcat.apache.org/tomcat-6.0-doc/index.html

2. Extract the FGB package, move the FGB package into 'webapps' folder 
or any other 'appBase' folder configured in server.xml of Tomcat.

3. Edit the config.xml.template to configure built-in tracks of the FGB. The 
config.xml file records each track's id, path, and type. 

Track IDs are the unique identifier to access the track. 

Track links should be the absolute path of track files. User can either 
specify one file to a track, or specify multiple files for different 
chromosomes of a track.

In the config.xml.template file, paths of attached annotations and databases
is '/fgbdata/....'. Users can simply put the attached 'fgbdata' folder in 
'/'. Then the absolute path of attached data in config.xml.template is no
need to be modified.

If the Track file has an index file, such as '.tbi', '.fai', etc., 
the index file must be put in the same folder as the track file, and 
the index file name must be 'full data file name' + standard suffix. Tomcat6 
must have the access to the track files, otherwise the tracks cannot 
be displayed in the FGB.

4. In addition, tracks in 'Basic' groups, including hg19 reference sequence,
cytoBand annotation, dbSNP database, HGNC genes, RefSeq Genes, Linkage
Disquilibriums data, and OMIM database are essential data/tool for the FGB, 
which are required to be prepared in prior. The FGB package data bundle includes
cytoBand, HGNC Genes, RefSeq Genes, Linkage Disquilibriums and OMIM data.
dbSNP142 and hg19 reference sequence must be downloaded before running the FGB.
	
a) dbSNP Build 142 (latest version)
ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi

dbSNP should be used in the 'Basic' track 'Snp'. In the config.xml.template 
we use a chr21 data from dbSNP_CEU population as demo.

b) hg19 Reference sequence
ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/
Please uncompress and merge all the Karyotype chromosome sequence into 
a single FASTA file. Then using software such as SAMtools 
(http://samtools.sourceforge.net/samtools.shtml#3) to generate a '.fai' 
index file.
"samtools faidx hg19.fa"

c) The FGB also hold a copy of above files for users' needs:
http://mlg.hit.edu.cn/FGB/data/dbsnp142.hg19.sorted.vcf.gz
http://mlg.hit.edu.cn/FGB/data/dbsnp142.hg19.sorted.vcf.gz.tbi
http://mlg.hit.edu.cn/FGB/data/hg19.tar.gz

Users can also integrate their own data into FGB by customizing the 
config.xml.template file.

5. After customizing the config.xml.template file, rename it to 'config.xml' 
and restart the Tomcat service to apply the changes.

6. Visit http://Your IP or Domain(:8080)/FGB/ to explore genomes in FGB.

Please feel free to contact pgbrowser@gmail.com for any further questions.

LJ Feb.3 2015
