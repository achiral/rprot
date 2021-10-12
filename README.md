# rbioc
- Docker RStudio sever with Bioconductor

# git clone
`git clone https://github.com/achiral/rbioc.git`

# change directory
`cd ./rbioc`

# edit template
## edit Dockerfile
- change Docker image (ex. `bioconductor/bioconductor_docker` => `FROM docker-ID/docker-image:latest`)  
- install R packages (see below exampe codes)  
`RUN install2.r pkgname1 pkgname2`  
`RUN Rscript -e "install.packages('vector_of_package_names')"`  
`RUN Rscript -e 'BiocManager::install("vector_of_package_names")'`  

## edit docker-compose.yml
- rename Docker image (ex. `image: docker-ID/docker-image:latest`)
## user approval
- change password password: `pw` => `<your_password>`
- delete environment option: `- DISABLE_AUTH=false` => `true` (without password sign in)

# docker compose
`docker-compose up --build -d`

# run RStudio from web browser
`localhost:8787`(`rstudio`/`pw` (default) or `<your_password>`)

# stop docker container
`docker-compose stop`

# restart docker container
`docker-compose start` without update
`docker-compose up --build -d` with update

# copy file from host to container
`docker cp <host directory> <container name>:<container ditectory>`

# copy file from container to host
`docker cp <container name>:<container directory> <host directory>`
