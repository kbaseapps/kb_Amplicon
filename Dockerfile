FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys FCAE2A0E115C3D8A
RUN echo 'deb https://cloud.r-project.org/bin/linux/debian stretch-cran35/' >> /etc/apt/sources.list

RUN apt-get update
RUN apt-get install -y r-base r-base-dev

RUN cp /usr/bin/R /kb/deployment/bin/.

## Install packages are available for ecologists
RUN Rscript -e "install.packages('vegan')"

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work

RUN mkdir -p /kb/module/work/vegan
RUN mkdir -p /kb/module/work/vegan/results

RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
