#!/bin/bash

git checkout --orphan newbranch
git add -A 
git commit -m 'cleanup'
git branch -D main
git branch -m main
git push -f origin main
git gc --aggressive --prune=all
