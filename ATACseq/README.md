
# Insert size scraper

```
for i in $(cat ../../samples.txt); do
SAMPLE=$i
MEDIAN_INSERT_SIZE=$(awk '/MEDIAN_INSERT_SIZE/{getline; print}' ${i}_insert_size_metrics.txt | cut -f 1,5)
echo -e "$SAMPLE\t$MEDIAN_INSERT_SIZE"
done > Insert_size_scraped.tsv
# add headder
echo -e 'sample\tMEDIAN_INSERT_SIZE\tMEAN_INSERT_SIZE' \
| cat - Insert_size_scraped.tsv > temp && mv temp Insert_size_scraped.tsv

lst Insert_size_scraped.tsv
```
