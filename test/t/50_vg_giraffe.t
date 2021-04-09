#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 20

vg construct -a -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg -G x.gbwt -v small/x.vcf.gz x.vg
vg snarls --include-trivial x.vg > x.snarls
vg index -s x.snarls -j x.dist x.vg
vg minimizer -k 29 -w 11 -g x.gbwt -i x.min x.xg

vg giraffe -x x.xg -H x.gbwt -m x.min -d x.dist -f reads/small.middle.ref.fq > mapped1.gam
is "${?}" "0" "a read can be mapped with all indexes specified without crashing"

vg minimizer -k 29 -b -s 18 -g x.gbwt -i x.sync x.xg

vg giraffe -x x.xg -H x.gbwt -m x.sync -d x.dist -f reads/small.middle.ref.fq > mapped.sync.gam
is "${?}" "0" "a read can be mapped with syncmer indexes without crashing"

rm -f x.vg x.xg x.gbwt x.snarls x.min x.sync x.dist x.gg

cp small/x.fa .
cp small/x.vcf.gz .
cp small/x.vcf.gz.tbi .

vg giraffe x.fa x.vcf.gz -f reads/small.middle.ref.fq > mapped2.gam
is "${?}" "0" "a read can be mapped with just FASTA and VCF without crashing"

# These files can differ as serialized and still represent the same data, due to protobuf field order not being specified.
# Tripping through JSON will sort all the keys.
vg view -aj mapped1.gam >mapped1.json
vg view -aj mapped2.gam >mapped2.json
vg view -aj mapped.sync.gam >mapped.sync.json

# Make sure at least one file converted successfully
SIZE="$(wc -c mapped2.json | cut -f1 -d' ')"
EMPTY=0
if [ "${SIZE}" == "0" ] ; then
    EMPTY=1
fi
is "${EMPTY}" "0" "mapping with just a FASTA and a VCF produced JSON-able alignments"

diff mapped1.json mapped2.json
is "${?}" "0" "mapping to manually-generated indexes and automatically-generated indexes is the same"

is "$(jq '.path' mapped1.json)" "$(jq '.path' mapped.sync.json)" "mapping with syncmers produces the same alignment as mapping with minimizers"

rm -rf mapped1.gam mapped1.json mapped2.gam mapped2.json mapped.sync.gam mapped.sync.json

vg giraffe x.fa x.vcf.gz -f small/x.fa_1.fastq > single.gam
is "$(vg view -aj single.gam | jq -c 'select((.fragment_next | not) and (.fragment_prev | not))' | wc -l)" "1000" "unpaired reads lack cross-references"

vg giraffe x.fa x.vcf.gz -f small/x.fa_1.fastq -f small/x.fa_1.fastq --fragment-mean 300 --fragment-stdev 100 > paired.gam
is "$(vg view -aj paired.gam | jq -c 'select((.fragment_next | not) and (.fragment_prev | not))' | wc -l)" "0" "paired reads have cross-references"

# Test paired surjected mapping
vg giraffe x.fa x.vcf.gz -iG <(vg view -a small/x-s13241-n1-p500-v300.gam | sed 's%_1%/1%' | sed 's%_2%/2%' | vg view -JaG - ) --output-format SAM >surjected.sam
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 4)" "$(printf '321\n762')" "surjection of paired reads to SAM yields correct positions"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 8)" "$(printf '762\n321')" "surjection of paired reads to SAM yields correct pair partner positions"
is "$(cat surjected.sam | grep -v '^@' | cut -f 1 | sort | uniq | wc -l)" "1" "surjection of paired reads to SAM yields properly matched QNAMEs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 7)" "$(printf '=\n=')" "surjection of paired reads to SAM produces correct pair partner contigs"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 2)" "$(printf '163\n83')" "surjection of paired reads to SAM produces correct flags"

# And unpaired surjected mapping
vg giraffe x.fa x.vcf.gz -G <(vg view -a small/x-s13241-n1-p500-v300.gam | sed 's%_1%/1%' | sed 's%_2%/2%' | vg view -JaG - ) --output-format SAM >surjected.sam
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 4)" "$(printf '321\n762')" "surjection of unpaired reads to SAM yields correct positions"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 8)" "$(printf '0\n0')" "surjection of unpaired reads to SAM yields correct pair partner positions"
is "$(cat surjected.sam | grep -v '^@' | cut -f 1 | sort | uniq | wc -l)" "2" "surjection of unpaired reads to SAM yields distinct QNAMEs"
is "$(cat surjected.sam | grep -v '^@' | cut -f 7)" "$(printf '*\n*')" "surjection of unpaired reads to SAM produces absent partner contigs"
is "$(cat surjected.sam | grep -v '^@' | sort -k4 | cut -f 2)" "$(printf '0\n16')" "surjection of unpaired reads to SAM produces correct flags"

rm -f x.vg x.gbwt x.gg x.snarls x.min x.dist x.gg x.fa x.fa.fai x.vcf.gz x.vcf.gz.tbi single.gam paired.gam surjected.sam

cp small/xy.fa .
cp small/xy.vcf.gz .
cp small/xy.vcf.gz.tbi .
vg giraffe xy.fa xy.vcf.gz -f small/x.fa_1.fastq -o SAM --ref-paths small/yx.dict | grep -E "^@(SQ|HD)" > surjected-yx.dict
vg giraffe xy.fa xy.vcf.gz -f small/x.fa_1.fastq -o SAM --ref-paths small/xy.dict | grep -E "^@(SQ|HD)" > surjected-xy.dict

diff surjected-yx.dict small/yx.dict
is "${?}" "0" "surjecting with a sequence dictionary in non-sorted order produces headers in non-sorted order"

diff surjected-xy.dict small/xy.dict
is "${?}" "0" "surjecting with a sequence dictionary in sorted order produces headers in sorted order"

rm -f xy.vg xy.gbwt xy.gg xy.snarls xy.min xy.dist xy.gg xy.fa xy.fa.fai xy.vcf.gz xy.vcf.gz.tbi surjected-yx.dict surjected-xy.dict

