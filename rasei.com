#!/bin/sh

fla=a3
flb=a4
num=32
rho=32

awk '$4 == "E" && $2%'$num' == 0 {print $1 / 1000.0 , $2/'$rho' }' sp.ras.$fla > sp.ras.$fla.xx.E
awk '$4 == "I" && $2%'$num' == 0 {print $1 / 1000.0 , $2/'$rho' }' sp.ras.$fla > sp.ras.$fla.xx.I

awk '$4 == "E" && $2%'$num' == 0 {print $1 / 1000.0 , $2/'$rho' }' sp.ras.$flb > sp.ras.$flb.xx.E
awk '$4 == "I" && $2%'$num' == 0 {print $1 / 1000.0 , $2/'$rho' }' sp.ras.$flb > sp.ras.$flb.xx.I

awk '{print $1 / 1000.0, $3}' sp.col.$fla > sp.col.$fla.xx.1
awk '{print $1 / 1000.0, $5}' sp.col.$fla > sp.col.$fla.xx.2

awk '{print $1 / 1000.0, $3}' sp.col.$flb > sp.col.$flb.xx.1
awk '{print $1 / 1000.0, $5}' sp.col.$flb > sp.col.$flb.xx.2

xmgrace -graph 0 sp.ras.$fla.xx.E \
        -graph 1 sp.ras.$flb.xx.E \
        -graph 2 sp.ras.$fla.xx.I \
        -graph 3 sp.ras.$flb.xx.I \
        -graph 4 sp.col.$fla.xx.1 \
        -graph 5 sp.col.$flb.xx.1 \
        -graph 6 sp.col.$fla.xx.2 \
        -graph 7 sp.col.$flb.xx.2 \
        -hdevice EPS -p rasei.gr  -printfile rasei.eps

/bin/rm sp.ras.$fla.xx.E sp.ras.$fla.xx.I
/bin/rm sp.ras.$flb.xx.E sp.ras.$flb.xx.I
/bin/rm sp.col.$fla.xx.1 sp.col.$fla.xx.2
/bin/rm sp.col.$flb.xx.1 sp.col.$flb.xx.2
