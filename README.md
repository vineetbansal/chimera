# chimera
A Web application + Python framework to find interaction sites for Protein Domains.

## Live Demo
This application is currently hosted for testing at https://protdomain.princeton.edu

### Installation

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

### Testing

The application can be run as a container using the podman/docker commands below.

```
podman build -t chimera .
```

Assuming `/path/to/pfam` is a folder with the structure outlined above, test the container by running the unit tests:

```
podman run -it --rm -v /path/to/pfam:/pfam:z localhost/chimera:latest pytest
```

The `-z` is an SELinux-specific flag. You may or may not need it.

### Running

If tests run okay:

Write a `.env` file at the root of the project that supplies envvars to the containers. At the time of this writing, it includes (values withheld):

```
SMTPHOST=
SMTPPORT=
SMTPUSER=
SMTPPASS=
PFAM_DIR=/path/to/pfam
```

`PFAM_DIR` should be set to whatever worked in the `Testing` section above.

Now run:

```
podman-compose up -d
```

The webapp runs on port `5000` but this can be controlled in `docker-compose.yml`.
