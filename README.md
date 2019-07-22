[![Build Status](https://travis-ci.org/vineetbansal/chimera.svg?branch=master)](https://travis-ci.org/vineetbansal/chimera)

# chimera
Web application + Python framework to find interaction sites for Protein Domains.

### Installation

```
docker build -t chimera .
docker run -it --rm chimera pytest
```

If tests run okay, start the Flask server (the default command).

```
docker run -it --rm -p 80:5000 chimera
```
