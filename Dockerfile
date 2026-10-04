# syntax=docker/dockerfile:1.7
# ROBIN application image: conda stack from robin.yml + pip-installed package.
# Nested tools (ClairS-To, MNP-Flex) still talk to a Docker daemon; mount the
# host socket and bind data at the same host paths. See docs/getting-started/docker.md.

FROM mambaorg/micromamba:2.3-debian12-slim

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        docker.io \
        git \
        git-lfs \
        procps \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER robin.yml /tmp/robin.yml
RUN micromamba install -y -n base -f /tmp/robin.yml \
    && micromamba clean --all --yes

WORKDIR /opt/robin
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml assets.json README.md ./
COPY --chown=$MAMBA_USER:$MAMBA_USER src ./src
COPY --chown=$MAMBA_USER:$MAMBA_USER docker ./docker

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG ROBIN_GIT_COMMIT=""
ENV ROBIN_GIT_COMMIT=${ROBIN_GIT_COMMIT} \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    XDG_CONFIG_HOME=/home/mambauser/.config

RUN chmod +x /opt/robin/docker/entrypoint.sh \
    && pip install --no-cache-dir .

# Prefer conda's libstdc++ so sqlite3/ICU (GUI auth) do not pick up Debian's older CXXABI.
ENV LD_LIBRARY_PATH=/opt/conda/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

# Pre-create so a named volume mounted here inherits mambauser ownership
# instead of being created as root (GUI security.db).
USER root
RUN mkdir -p /home/$MAMBA_USER/.config/robin /home/$MAMBA_USER/.local/robin \
    && chown -R $MAMBA_USER:$MAMBA_USER /home/$MAMBA_USER/.config /home/$MAMBA_USER/.local
USER $MAMBA_USER

EXPOSE 8081 8265

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/opt/robin/docker/entrypoint.sh"]
CMD ["robin", "--help"]
