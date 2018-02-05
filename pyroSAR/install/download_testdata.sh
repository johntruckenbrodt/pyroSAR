#!/usr/bin/env bash


mkdir -p $TESTDATA_DIR

#cd $TESTDATA_DIR

echo "Start Download forest_brazil"
wget --quiet -P $TESTDATA_DIR 'ftp://ftp.eorc.jaxa.jp/pub/ALOS-2/1501sample/310_forestbrazil/0000022708_001001_ALOS2015976960-140909.zip'
echo "End download forest_brazil"
