#!/usr/bin/env bash
cd $H

mkdir Filter

cat ./preliminary_ABC/proteomes/* >> ./Filter/ABC_preliminary_total.faa
/home/pioannidis/Programs/interproscan-5.30-69.0/interproscan.sh -appl pfam -i ./Filter/ABC_preliminary_total.faa -o ./Filter/IPSCAN.tsv -f TSV


