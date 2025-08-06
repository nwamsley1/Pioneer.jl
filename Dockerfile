FROM ubuntu:24.04

RUN apt-get update && apt-get install -y curl ca-certificates jq && rm -rf /var/lib/apt/lists/*

# Get latest release info from GitHub API and download the correct .deb asset
RUN LATEST_VERSION=$(curl -s https://api.github.com/repos/nwamsley1/Pioneer.jl/releases/latest | jq -r .tag_name) && \
    FILE_NAME=$(curl -s https://api.github.com/repos/nwamsley1/Pioneer.jl/releases/latest \
        | jq -r '.assets[] | select(.name | endswith(".deb")) | .name') && \
    curl -L "https://github.com/nwamsley1/Pioneer.jl/releases/download/${LATEST_VERSION}/${FILE_NAME}" -o /tmp/${FILE_NAME} && \
    apt-get install -y /tmp/${FILE_NAME} && \
    rm /tmp/${FILE_NAME}

ENV PATH="/usr/local/pioneer/bin:$PATH"

RUN pioneer --help
