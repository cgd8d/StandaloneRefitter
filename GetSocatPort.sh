#!/bin/bash

# WARNING: this script has a race condition.  Do not let multiple scripts call it near together in time.

USEPORT=13967
while [ 1 ]; do
  INUSE=`python -c "import socket; s = socket.socket(); s.bind(('0.0.0.0', $USEPORT))" > /dev/null 2>&1 && echo "open" || echo "inuse"`
  if [ x"$INUSE" = "xopen" ]; then
    break
  fi
  USEPORT=$(($USEPORT + 1))
done
echo $USEPORT
