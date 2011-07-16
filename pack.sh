#!/bin/bash

bash build.sh
tar cvfz team5-submission.tar.gz \
    pack.sh run.sh build.sh \
    team5.jar team5.jardesc team5.jar.manifest .classpath .project \
    bin \
    src \
    cpp \
    bin-cpp
    
