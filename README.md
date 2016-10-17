###gatk-tools: Utilities for processing sequencing data and genomic variants using GATK

This package contains some utilities for processing sequencing data and genomic
variants using the framework provided by the
[Genome Analysis Toolkit](http://www.broadinstitute.org/gatk) (GATK). 

####Installing and running gatk-tools

Please see [README](docs/README) file for more information on the utilities available
within the package and details of how to install and run the tools from a
[pre-packaged release](https://github.com/crukci-bioinformatics/gatk-tools/releases).

For more information on the utilities available in this package, and how to
install and run these tools, please see the [README](docs/README) file.

####Building from source

The gatk-tools package is built using Apache Maven, a software project
management and build automation tool. Details on how to install and run Maven
can be found [here](http://maven.apache.org).

Maven will automatically download dependencies from the central Maven repository
with the exception of the GATK library. GATK is available under a dual licensing
model with separate licenses for academic, non-commerical research and
for-profit purposes, and can be downloaded from the
[GATK site](http://www.broadinstitute.org/gatk). The GATK JAR file will need to
be installed separately as shown in the following instructions.

1. Clone the project

        git clone https://github.com/crukci-bioinformatics/gatk-tools.git
        cd gatk-tools

2. Install the GATK jar file to the local Maven repository (substituting the
version and path to the GATK jar as appropriate)

        mvn install:install-file -Dpackaging=jar -DgroupId=org.broadinstitute.gatk -DartifactId=gatk -Dversion=3.6 -Dfile=GenomeAnalysisTK.jar

3. Build and package the gatk-tools package

        mvn package

4. Unpack the gatk-tools tarball to an installation directory (substituting the
version number as appropriate)

        tar zxf target/gatk-tools-1.0-distribution.tar.gz

This will create a directory named gatk-tools-1.0 which can be moved to the
desired installation location.

