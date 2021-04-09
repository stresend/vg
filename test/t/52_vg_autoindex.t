#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 21

vg autoindex -p auto -w map -r tiny/tiny.fa -v tiny/tiny.vcf.gz --force-unphased
is $(echo $?) 0 "autoindexing successfully completes indexing for vg map with basic input"
is $(ls auto.* | wc -l) 3 "autoindexing makes 3 outputs for vg map" 
is $(ls auto.xg | wc -l) 1 "autoindexing makes an XG for vg map"
is $(ls auto.gcsa* | wc -l) 2 "autoindexing makes a GCSA2/LCP pair for vg map"
vg sim -x auto.xg -n 20 -a -l 10 | vg map -d auto -t 1 -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by vg map"

rm auto.*

vg autoindex -p auto -w map -r small/x.fa -v small/x.vcf.gz -r small/y.fa -v small/y.vcf.gz 
is $(echo $?) 0 "autoindexing successfully completes indexing for vg map with chunked input"
vg sim -x auto.xg -n 20 -a -l 10 | vg map -d auto -t 1 -G - > /dev/null
is $(echo $?) 0 "chunked autoindexing results can be used by vg map"

rm auto.*

vg autoindex -p auto -w map -r small/xy.fa -v small/xy2.vcf.gz
is $(echo $?) 0 "autoindexing successfully completes indexing for vg map with phased input"
vg sim -x auto.xg -n 20 -a -l 10 | vg map -d auto -t 1 -G - > /dev/null
is $(echo $?) 0 "phased autoindexing results can be used by vg map"

rm auto.*

# to add: GFA construction

vg autoindex -p auto -w mpmap -r tiny/tiny.fa -v tiny/tiny.vcf.gz -x tiny/tiny.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with unchunked input"
is $(ls auto.* | wc -l) 6 "autoindexing creates 6 files for mpmap/rpvg"
vg sim -x auto.spliced.xg -n 20 -a -l 10 | vg mpmap -x auto.spliced.xg -g auto.spliced.gcsa -d auto.spliced.dist -B -t 1 -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by vg mpmap"
is $(vg paths -g auto.haplotx.gbwt -L | wc -l) 4 "haplotype transcript GBWT made by autoindex is valid"
is $(cat auto.txorigin.tsv | wc -l) 5 "transcript origin table has expected number of rows" 

rm auto.*

vg autoindex -p auto -w mpmap -r small/x.fa -r small/y.fa -v small/x.vcf.gz -v small/y.vcf.gz -x small/x.gtf -x small/y.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with chunked input"
is $(ls auto.* | wc -l) 6 "autoindexing creates 6 files for mpmap/rpvg with chunked input"

rm auto.*

vg autoindex -p auto -w mpmap -r small/x.fa -r small/y.fa -v small/x.vcf.gz -v small/y.vcf.gz -x small/xy.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with partially chunked input"

rm auto.*

vg autoindex -p auto -w mpmap -r small/xy.fa -v small/x.vcf.gz -v small/y.vcf.gz -x small/xy.gtf
is $(echo $?) 0 "autoindexing successfully completes indexing for vg mpmap with another partially chunked input"

rm auto.*

vg autoindex -p auto -w giraffe -r tiny/tiny.fa -v tiny/tiny.vcf.gz 
is $(echo $?) 0 "autoindexing successfully completes indexing for vg giraffe with unchunked input"
is $(ls auto.* | wc -l) 4 "autoindexing creates 4 inputs for vg giraffe"
vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz > t.vg
vg index -x t.xg t.vg
vg sim -x t.xg -n 20 -a -l 10 | vg giraffe -g auto.gg -H auto.giraffe.gbwt -m auto.min -d auto.dist -G - > /dev/null
is $(echo $?) 0 "basic autoindexing results can be used by vg giraffe"

rm auto.*
rm t.*



