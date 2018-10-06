UPDATE: The Family Genome Browser service has been moved to http://www.pgbrowser.org/FGB/

The installation package can be downloaded at http://www.pgbrowser.org/FGB/FGB-1.0.tar.gz

The FGB data bundle can be downloaded at http://www.pgbrowser.org/FGB/FGB_data_bundle.tar.gz

LJ Oct.6 2018

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
which are required to be prepared in prior. The FGB package data bundle is
available at http://mlg.hit.edu.cn/FGB/FGB_data_bundle.tar.gz
	
Users can also integrate their own data into FGB by customizing the 
config.xml.template file.

5. After customizing the config.xml.template file, rename it to 'config.xml' 
and restart the Tomcat service to apply the changes.

6. Visit http://Your IP or Domain(:8080)/FGB/ to explore genomes in FGB.

Please feel free to contact pgbrowser@gmail.com for any further questions.

LJ Feb.13 2015
