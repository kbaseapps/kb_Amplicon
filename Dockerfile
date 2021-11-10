FROM kbase/sdkbase2:python
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# R related installations
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys FCAE2A0E115C3D8A
RUN echo 'deb https://cloud.r-project.org/bin/linux/debian stretch-cran35/' >> /etc/apt/sources.list

RUN apt-get update --fix-missing
RUN apt-get install -y r-base r-base-dev

RUN cp /usr/bin/R /kb/deployment/bin/.
RUN cp /usr/bin/Rscript /kb/deployment/bin/.

## Install packages are available for ecologists
# vegan: Community Ecology Package
RUN Rscript -e "install.packages('vegan')"
# vegan3d: Static and Dynamic 3D Plots for the 'vegan' Package
RUN Rscript -e "install.packages('vegan3d')"

## Install other packages
# parmigene: Parallel Mutual Information estimation for Gene Network reconstruction
RUN Rscript -e "install.packages('parmigene')"

# jsonlite: A Robust, High Performance JSON Parser and Generator for R
RUN Rscript -e "install.packages('jsonlite')"

RUN pip install --upgrade pip \
    && python --version

RUN pip install coverage==6.1.1 \
    && pip install pandas==1.1.5 \
    && pip install dotmap==1.3.25 \
    && pip install plotly==5.3.1 \
    && pip install mock==4.0.3

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work

RUN mkdir -p /kb/module/work/amplicon
RUN mkdir -p /kb/module/work/amplicon/results

RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
