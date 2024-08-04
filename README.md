# webapps

These are several Pluto notebooks bundled to run as webapps through a container. The apps can be viewed live here via the University of California, Riverside notebook server: [https://notebooks.engr.ucr.edu](https://notebooks.engr.ucr.edu).   


## Build container 

1. Clone via git

```bash
$ git clone https://github.com/mdpetters/webapps.git
```

2. Build container

```bash
$ cd webapps/virtualDMA
$ docker build . -t docker.io/mdpetters/virtualtdma:latest
```

You can replace `docker` with `podman`. You can replace the prefix `docker.io/mdpetters/` as desired.

