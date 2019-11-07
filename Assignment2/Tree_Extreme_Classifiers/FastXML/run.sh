cat train.X | sed -e 's/^[ \t]*//' > temp
mv temp train.X
cat train.y | sed -e 's/^[ \t]*//' > temp
mv temp train.y
cat test.X | sed -e 's/^[ \t]*//' > temp
mv temp test.X
cat test.y | sed -e 's/^[ \t]*//' > temp
mv temp test.y