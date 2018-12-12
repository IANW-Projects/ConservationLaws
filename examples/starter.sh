#!/bin/sh

cat ../python/*.py > runfile.py
cat ../python/BalanceLaws/*.py >> runfile.py
cp ../python/cl/cl.py cl.py
cat $1 >> runfile.py
python3 runfile.py
rm cl.py
rm runfile.py

