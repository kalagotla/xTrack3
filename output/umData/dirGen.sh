n=177;
max=573;
while [ "$n" -le "$max" ]; do
  mkdir "$n"
  n=`expr "$n" + 1`;
done
