if [ "$#" -ne 4 ]; then
    echo "usage: analysis.sh nb_dir tree_dir depth result_dir $#"
    exit
fi
nb_dir=$1
tree_dir=$2
depth=$3
result_dir=$4
cat "" > $nb_dir/info_raw
cat "" > $nb_dir/info
for i in $(seq 0 $depth); do echo $i; echo $i >> $nb_dir/info_raw; grep -o ' ' $nb_dir/$i  | wc -l >> $nb_dir/info_raw; done
cat $nb_dir/info_raw | paste - - > $nb_dir/info
cat $nb_dir/info | awk '{print $2}' > $result_dir/nb_count_raw
cat $tree_dir/info | awk '{print $2}' > $result_dir/tree_count
paste $result_dir/nb_count_raw $result_dir/tree_count | head -n $depth | awk '{print $1 - $2}' > $result_dir/nb_count
paste $result_dir/nb_count $result_dir/tree_count | head -n $depth | awk '{print $1 / $2}' > $result_dir/div
