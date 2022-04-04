ARG BASE_CONTAINER=sagemath/sagemath:9.3
FROM $BASE_CONTAINER

# USER root
# COPY ./entrypoint-dev.sh /usr/local/bin/sage-entrypoint

ARG HOME=/home/sage
USER sage
ENV HOME $HOME
WORKDIR $HOME

ENTRYPOINT [ "/usr/local/bin/sage-entrypoint" ]
CMD ["bash"]
