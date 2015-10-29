
To build the Docker image, move to the folder with the desired configuration and run Docker build:

    cd cuda_331_openblas/
    docker build --no-cache -t gadgetron .

The `--no-cache` option will ensure that you do a fresh pull from the git repos. To start Docker container from image:

    export CUDA_DEVICES="--device=/dev/nvidia0:/dev/nvidia0 --device=/dev/nvidiactl:/dev/nvidiactl --device=/dev/nvidia-uvm:/dev/nvidia-uvm"
    docker run -e "GADGETRON_RELAY_HOST=<MY RELAY HOST IP/NAME>"  ${CUDA_DEVICES} --name gt1 -p 9002:9002 --rm -t gadgetron

You will now have a container named `gt1` running and the Gadgetron inside will be exposed on port 9002 on the Docker host. The `CUDA_DEVICES` part is only needed if you would like to expose your CUDA devices inside the Docker container. There is a helper script `nvidia_devices.sh`, which will generate this string for you. 
