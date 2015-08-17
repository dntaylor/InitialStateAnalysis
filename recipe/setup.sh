#!/usr/bin/env bash

pushd $CMSSW_BASE/src

git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
pushd HiggsAnalysis/CombinedLimit
git checkout 74x-root6
popd

popd
