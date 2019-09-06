[![Build Status](https://travis-ci.org/vineetbansal/chimera.svg?branch=master)](https://travis-ci.org/vineetbansal/chimera)

# chimera
A Web application + Python framework to find interaction sites for Protein Domains.

## Live Demo
This application is currently hosted for testing at http://protdomain.com

Note that during development, server resources allocated to the above website are minimal, since this is only for development purposes. Queries will take substantially longer than they would in the final release. 

### Installation

The application can be run as a container using the Docker commands below.

```
docker build -t chimera .
docker run -it --rm chimera pytest
```

If tests run okay, start the Flask server (the default command).
```
docker run -d --rm -p 80:5000 -v /path/to/pfam/files:/pfam chimera
```

Note that a valid [Pfam Database](https://pfam.xfam.org/) needs to be accessible *from within* the container. Pfam files are not included in this source (or as part of the built container image) because of size considerations.

Specifically, `Pfam-A.hmm.dat` and `Pfam-A.hmm` files (with `hmmpress` run on the latter) are needed at the location `/pfam` within the container.

This can be accomplished by downloading the required files from [ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release) on the host machine, running `hmmpress` on it, and then mounting the folder a volume within the container.

On the demo server above, the contents of the `/mnt/pfam_data` look like the following:

```
-rw-r--r-- 1 root root 1434541235 Aug 23 14:31 Pfam-A.hmm
-rw-r--r-- 1 root root    2935721 Aug 23 14:31 Pfam-A.hmm.dat
-rw-r--r-- 1 root root  328550912 Aug 23 14:33 Pfam-A.hmm.h3f
-rw-r--r-- 1 root root    1237206 Aug 23 14:33 Pfam-A.hmm.h3i
-rw-r--r-- 1 root root  593831723 Aug 23 14:33 Pfam-A.hmm.h3m
-rw-r--r-- 1 root root  698539111 Aug 23 14:33 Pfam-A.hmm.h3p
```
The container is then run with the options `-v /mnt/protdomain_data:/pfam` to make this folder accessible to the container.
