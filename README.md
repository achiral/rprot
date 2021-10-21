# rbioc
- Docker RStudio sever with Bioconductor

# git clone
`git clone https://github.com/achiral/rprot.git`

# change directory
`cd ./rprot`

# edit template
## edit Dockerfile
- change Docker image (ex. `achiral/rbioc:latest` => `FROM docker-ID/docker-image:latest`)  
- install R packages (see below exampe codes)  
`RUN install2.r pkgname1 pkgname2`
`RUN Rscript -e "install.packages('vector_of_package_names')"`
`RUN Rscript -e 'BiocManager::install("vector_of_package_names")'`

## edit docker-compose.yml
- rename Docker image (ex. `image: docker-ID/docker-image:latest`)
- rename Docker container (ex. `container_name: docker-ID/docker-container:latest`)
## user approval
- change password (ex. `pw` => `<your_password>`)
- change environment option: (`- DISABLE_AUTH=true` => `false` (activate password sign in))

# docker compose
`docker-compose up --build -d`

# run RStudio from web browser
`open http://localhost:8787/`(`rstudio`/`pw` (default) or `<your_password>`)

# stop docker container
`docker-compose stop` (or `docker-compose down`)

# restart docker container
`docker-compose start` without update (stopped docker container with `docker-compose stop`)
`docker-compose up --build -d` with update (stopped docker container with `docker-compose down`)

# copy file from host to container
`docker cp <host directory> <container name>:<container ditectory>`
(ex. `docker cp /Users/user/Dropbox/data/. rprot:/home/rstudio/rproject/data/`)

# copy file from container to host
`docker cp <container name>:<container directory> <host directory>`
