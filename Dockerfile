ARG BASE_CONTAINER=sagemath/sagemath:9.3
FROM $BASE_CONTAINER

USER root
ENTRYPOINT [ "/usr/local/bin/sage-entrypoint" ]
