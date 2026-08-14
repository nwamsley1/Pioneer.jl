FROM ubuntu:24.04

# The apt lists are deliberately kept until the .deb is installed below. They
# used to be deleted here, which left the next layer with no package index: the
# .deb's dependencies then reported as "not installable" even though they are in
# noble's main component. That was latent until the package first declared any
# Depends at all, and it failed the v0.8.0 Docker publish.
RUN apt-get update \
    && apt-get install -y curl ca-certificates jq libicu-dev

ARG PIONEER_VERSION=latest

# Get release info from GitHub API and download the correct .deb asset
RUN if [ "$PIONEER_VERSION" = "latest" ]; then \
        RELEASE_URL="https://api.github.com/repos/nwamsley1/Pioneer.jl/releases/latest"; \
    else \
        RELEASE_URL="https://api.github.com/repos/nwamsley1/Pioneer.jl/releases/tags/${PIONEER_VERSION}"; \
    fi && \
    TAG_NAME=$(curl -s "$RELEASE_URL" | jq -r .tag_name) && \
    FILE_NAME=$(curl -s "$RELEASE_URL" \
        | jq -r '.assets[] | select(.name | endswith(".deb")) | .name') && \
    curl -L "https://github.com/nwamsley1/Pioneer.jl/releases/download/${TAG_NAME}/${FILE_NAME}" -o /tmp/${FILE_NAME} && \
    apt-get update && \
    apt-get install -y /tmp/${FILE_NAME} && \
    rm /tmp/${FILE_NAME} && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="/usr/local/pioneer/bin:$PATH"

RUN pioneer search --help
