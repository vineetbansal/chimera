[![Build Status](https://travis-ci.org/vineetbansal/chimera.svg?branch=master)](https://travis-ci.org/vineetbansal/chimera)

# chimera
Minimal Flask+Docker Boilerplate

### Installation

```
docker build -t chimera .
docker run -it --rm chimera pytest
```

If tests run okay, start the Flask server (the default command).

```
docker run -it --rm -p 80:5000 chimera
```
