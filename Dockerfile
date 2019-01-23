# use a base image that has gcc and cmake and fix the version
# so that everyone has the same
FROM spritsail/alpine-cmake:3.6

WORKDIR /app

ADD . /app

# generate the makefile, compile the project and run unit tests
RUN cmake -G "Unix Makefiles" && \
    cmake --build . && \
    cd tests && \
    ctest && \
    cd ../

# start the TestU01 tests when the container is run
CMD [ "./prog" ]


# Here we use docker in order everyone to have the
# same execution environement to reproduce our tests
# by running our code.