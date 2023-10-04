FROM rocker/rstudio

RUN /rocker_scripts/install_tidyverse.sh

ENV CTAN_REPO=https://mirror.ctan.org/systems/texlive/tlnet
ENV PATH=$PATH:/usr/local/texlive/bin/linux

RUN /rocker_scripts/install_verse.sh

COPY ./install.sh /home/rstudio
RUN chmod +x /home/rstudio/install.sh
RUN /home/rstudio/install.sh

COPY . /home/rstudio

EXPOSE 8787

CMD [ "/init" ]