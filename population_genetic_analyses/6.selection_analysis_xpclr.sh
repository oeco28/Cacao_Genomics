# . for each chromosome

for((i=2;i<=10;i++)); do
        ./XPCLR -xpclr criollo.chr$i curaray.chr$i chr$i.map cacao.chr$i.xpclr.out -w1 0.005 200 1000 1 -p0 0.95
;
done
