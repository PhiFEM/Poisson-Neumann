#!/bin/bash
${CONTAINER_ENGINE} run --privileged -ti -v $(pwd):/data -w /data docker.io/freefem/freefem:latest bash