FROM alpine:3.14
RUN apk add boost-dev
RUN apk add qt5-qtbase qt5-qtbase-dev
RUN apk add --update npm
RUN apk add nodejs
RUN apk add git
RUN apk add g++
RUN apk add gcc
RUN apk add make
RUN apk add cmake
RUN apk add eigen-dev

COPY build_unsqlite.sh /home/ubuntu/
RUN ls /home/ubuntu/
RUN /home/ubuntu/build_unsqlite.sh


# --no-cache

